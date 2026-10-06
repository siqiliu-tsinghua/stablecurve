// 微基准：GWloop vs fastGWloop（m=6，全部同 idx=1 组）
fn main() {
    use sc_gwall::gwloop::{LoopIns, fast_gwloop, gwloop};
    // 场景 2：跨 3 组（idx 1,1,2,2,3,3）——fast 枚举 (2!)³ vs slow 的 6!
    let ta: Vec<LoopIns> = vec![
        LoopIns { tau: 0, idx: 1 },
        LoopIns { tau: 1, idx: 1 },
        LoopIns { tau: 0, idx: 2 },
        LoopIns { tau: 1, idx: 2 },
        LoopIns { tau: 0, idx: 3 },
        LoopIns { tau: 1, idx: 3 },
    ];
    let tb: Vec<LoopIns> = vec![
        LoopIns { tau: 1, idx: 1 },
        LoopIns { tau: 0, idx: 1 },
        LoopIns { tau: 1, idx: 2 },
        LoopIns { tau: 0, idx: 2 },
        LoopIns { tau: 1, idx: 3 },
        LoopIns { tau: 0, idx: 3 },
    ];
    let t0 = std::time::Instant::now();
    let a = gwloop(3, 2, &[0, 0], &ta, &tb);
    let t1 = t0.elapsed();
    let t2 = std::time::Instant::now();
    let b = fast_gwloop(3, 2, &[0, 0], &ta, &tb);
    let t3 = t2.elapsed();
    println!("gwloop: {} ({t1:?})", a);
    println!("fast  : {} ({t3:?})", b);
    assert_eq!(a.to_string(), b.to_string());
}
