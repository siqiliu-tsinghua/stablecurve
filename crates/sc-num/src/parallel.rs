//! 有界的作用域并行：调用线程参与计算，嵌套调用自动串行。
//! 作业按固定块编号领取，结果按块序收集；工作状态由各线程独占。

use std::cell::Cell;
use std::sync::atomic::{AtomicUsize, Ordering};

thread_local! {
    static ACTIVE: Cell<bool> = const { Cell::new(false) };
}

pub fn thread_count(variable: &str, default_max: usize) -> usize {
    if ACTIVE.with(Cell::get) {
        return 1;
    }
    let available = std::thread::available_parallelism().map_or(1, usize::from);
    std::env::var(variable)
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .filter(|&n| n > 0)
        .unwrap_or(default_max)
        .min(available)
        .max(1)
}

/// 在外层作业中禁止再次创建线程，退出时恢复原状态（包括 panic 展开）。
pub fn serial_scope<T>(f: impl FnOnce() -> T) -> T {
    struct Reset(bool);
    impl Drop for Reset {
        fn drop(&mut self) {
            ACTIVE.with(|v| v.set(self.0));
        }
    }
    let _reset = Reset(ACTIVE.with(|v| v.replace(true)));
    f()
}

/// 每个工作线程初始化一次局部状态，按需领取固定编号的块。
/// count=0 不创建线程；threads=1 也使用相同的块边界和结果顺序。
pub fn chunks<S, R: Send>(
    count: usize,
    chunk: usize,
    threads: usize,
    init: impl Fn() -> S + Sync,
    work: impl Fn(&mut S, std::ops::Range<usize>) -> R + Sync,
) -> Vec<R> {
    assert!(chunk > 0);
    let blocks = count.div_ceil(chunk);
    if blocks == 0 {
        return Vec::new();
    }
    let threads = if ACTIVE.with(Cell::get) {
        1
    } else {
        threads.max(1).min(blocks)
    };
    if threads == 1 {
        return serial_scope(|| {
            let mut state = init();
            (0..blocks)
                .map(|b| work(&mut state, b * chunk..((b + 1) * chunk).min(count)))
                .collect()
        });
    }
    let next = AtomicUsize::new(0);
    let run = || {
        serial_scope(|| {
            let mut state = init();
            let mut results = Vec::new();
            loop {
                let b = next.fetch_add(1, Ordering::Relaxed);
                if b >= blocks {
                    break;
                }
                results.push((b, work(&mut state, b * chunk..((b + 1) * chunk).min(count))));
            }
            results
        })
    };
    let mut results = std::thread::scope(|scope| {
        let workers: Vec<_> = (1..threads).map(|_| scope.spawn(run)).collect();
        let mut results = run();
        for worker in workers {
            match worker.join() {
                Ok(part) => results.extend(part),
                Err(payload) => std::panic::resume_unwind(payload),
            }
        }
        results
    });
    results.sort_unstable_by_key(|(b, _)| *b);
    results.into_iter().map(|(_, value)| value).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_item_once_in_order_and_nested_work_stays_serial() {
        for threads in [1, 2, 4] {
            let visits: Vec<_> = (0..1001).map(|_| AtomicUsize::new(0)).collect();
            let values = chunks(
                visits.len(),
                13,
                threads,
                || (),
                |_, range| {
                    assert_eq!(thread_count("SC_NUM_TEST_UNSET_THREADS", 16), 1);
                    range
                        .map(|i| {
                            visits[i].fetch_add(1, Ordering::Relaxed);
                            i * i
                        })
                        .collect::<Vec<_>>()
                },
            )
            .into_iter()
            .flatten()
            .collect::<Vec<_>>();
            assert_eq!(values, (0..visits.len()).map(|i| i * i).collect::<Vec<_>>());
            assert!(visits.iter().all(|v| v.load(Ordering::Relaxed) == 1));
        }
    }
}
