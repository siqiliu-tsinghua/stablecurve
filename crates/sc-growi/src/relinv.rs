//! 相对理论的多分量分裂枚举器——原版 `libgw/relinv.{h,cc}` 的逐行直译。
//!
//! 这是全程序最复杂的状态机：枚举度/亏格/类在任意多个分量上的全部分配、
//! 环路节点数、节点类与节点重数，并自动去除对称构型（same 值块跟踪，
//! 对称因子折进 [`RelInv::cof`]）。数据布局与原版完全一致：
//! - `value[m][n]`：分量 m 上第 n 类的个数；n=0 度、n=1 亏格+节点−1、
//!   n=2 环路节点数；n≥3 为 distinct 类（类列表 `cls`）
//! - `same[m][n]`：只考虑前 n 层数据时，从 m 开始有多少个相邻分量相同
//! - `node_{class,mult,same}[m][j]`：分量 m 第 j 个节点的类/重数/节点同值
//!
//! M1 只用亏格 0（relsplit），节点平面按 genus+1 分配，亏格 0 即单节点。

use crate::inv::{Inv, Pdata};
use crate::num::{One, Q, Zero, q};

use crate::split2::factorial;

// —— 下标访问宏（直译原版的 value()/same()/node_*() 内联函数）——
// 索引用 i32 保存（与原版 int 一致），访问处转换。

macro_rules! value {
    ($s:ident, $m:expr, $n:expr) => {
        $s.value[($n as usize) * $s.firstmax + ($m as usize)]
    };
}
macro_rules! same {
    ($s:ident, $m:expr, $n:expr) => {
        $s.same[($n as usize) * $s.firstmax + ($m as usize)]
    };
}
macro_rules! node_mult {
    ($s:ident, $m:expr, $n:expr) => {
        $s.node[(3 * ($n as usize)) * $s.firstmax + ($m as usize)]
    };
}
macro_rules! node_class {
    ($s:ident, $m:expr, $n:expr) => {
        $s.node[(3 * ($n as usize) + 1) * $s.firstmax + ($m as usize)]
    };
}
macro_rules! node_same {
    ($s:ident, $m:expr, $n:expr) => {
        $s.node[(3 * ($n as usize) + 2) * $s.firstmax + ($m as usize)]
    };
}
macro_rules! cdeg_ {
    ($s:ident, $m:expr) => {
        $s.value[($m as usize)]
    };
}
macro_rules! cgenus_ {
    ($s:ident, $m:expr) => {
        $s.value[$s.firstmax + ($m as usize)]
    };
}
macro_rules! cloop_ {
    ($s:ident, $m:expr) => {
        $s.value[2 * $s.firstmax + ($m as usize)]
    };
}

pub struct RelInv {
    // 从 iv 拷贝的常量
    #[allow(dead_code)]
    points: i32,
    totalgenus: i32,
    xdim: i32,
    xdeg: i32,

    // distinct 类（前三个槽位保留：0=度 1=亏格+节点 2=环路）
    ncls: i32,
    cls: Vec<Pdata>,

    ncpt: i32,
    skip: bool,
    skipclass: Pdata,
    noclass: bool,

    /// 原始类的个数（欧拉示性数公式；有理数因为中间量是分数幂）
    prims: Q,
    /// 类分裂的组合因子（多项式系数），随分裂变化持续更新
    coeff: Q,
    node_minclass: i32,

    cls_vdim: Vec<i32>,
    cpt_vdim: Vec<i32>,
    cpt_tm: Vec<i32>,
    cpt_prim: Vec<i32>,

    value: Vec<i32>,
    same: Vec<i32>,
    firstmax: usize,
    node: Vec<i32>,
}

impl RelInv {
    /// 构造。`skipc` = 必须放在内部分量上的类下标（-1 = 无）。
    pub fn new(iv: &Inv, skipc: i32) -> Self {
        let firstmax = (iv.cdeg + 1) as usize;

        let mut r = RelInv {
            points: iv.pts.len() as i32,
            totalgenus: iv.genus as i32,
            xdim: iv.dim as i32,
            xdeg: iv.deg as i32,
            ncls: 3,
            cls: vec![Pdata::default(); iv.pts.len() + 3],
            ncpt: 0,
            skip: false,
            skipclass: Pdata::default(),
            noclass: false,
            prims: Q::ZERO,
            coeff: Q::ONE,
            node_minclass: 0,
            cls_vdim: vec![0; iv.pts.len() + 3],
            cpt_vdim: vec![0; firstmax],
            cpt_tm: vec![0; firstmax],
            cpt_prim: vec![0; firstmax],
            value: vec![0; (iv.pts.len() + 6) * firstmax],
            same: vec![0; (iv.pts.len() + 8) * firstmax],
            firstmax,
            // 原版按 3*(genus+1) 分配，但 nodeclass_next 在 ni=cloop 时会写
            // ni+1 平面（g=0 时被守卫挡住）；多留一层防御
            node: vec![0; 3 * (iv.genus as i32 + 2).max(2) as usize * firstmax],
        };

        // 读入 distinct 类
        for i in 0..iv.pts.len() as i32 {
            if i == skipc {
                r.skip = true;
                r.skipclass = iv.pts[i as usize];
            } else if i == 0 || iv.pts[i as usize] != iv.pts[i as usize - 1] {
                let ncls = r.ncls as usize;
                r.cls[ncls] = iv.pts[i as usize];
                r.cls_vdim[ncls] = -iv.pts[i as usize].vdim(iv.dim);
                value!(r, 0, r.ncls) = 1;
                if (iv.pts[i as usize].h as i32) < r.cpt_prim[0] {
                    r.cpt_prim[0] = iv.pts[i as usize].h as i32;
                }
                r.ncls += 1;
            } else {
                value!(r, 0, r.ncls - 1) += 1;
            }
        }

        // 初始：只有内部分量，最大度与亏格，全部类在内
        r.ncpt = 1;
        cdeg_!(r, 0) = iv.cdeg as i32;
        cgenus_!(r, 0) = iv.genus as i32;
        cloop_!(r, 0) = 0;

        // 亏格变化 1 的虚维数代价
        r.cls_vdim[1] = 2 * (2 - iv.dim as i32);

        for i in 0..r.ncls {
            same!(r, 0, i) = 1;
        }

        // 原始类个数：欧拉示性数 = top 陈类 = (N+1) + (-1)^N · #primitive
        if iv.deg != 1 {
            let mut p = (q_pow(q(1 - iv.deg as i64), iv.dim as i64 + 2) + q(iv.deg as i64 - 1))
                / q(iv.deg as i64);
            if iv.dim & 1 != 0 {
                p = -p;
            }
            r.prims = p;
        }

        // 节点最小 H 类（g=0 无点时的加速剪枝）
        r.node_minclass = 0;
        if r.totalgenus == 0 && iv.pts.is_empty() {
            if r.xdeg == r.xdim + 2 {
                r.node_minclass = r.xdim - 1; // Calabi-Yau
            }
            if r.xdeg == r.xdim + 1 {
                r.node_minclass = r.xdim; // Fano 指标 1
            }
            if r.xdeg < r.xdim + 1 {
                r.node_minclass = r.xdim + 1; // Fano 指标 >1
            }
        }

        r.noclass = iv.pts.is_empty() || (iv.pts.len() == 1 && r.skip);

        r
    }

    /// 外部分量个数。
    #[inline]
    pub fn cpts(&self) -> usize {
        (self.ncpt - 1) as usize
    }

    // —— deg_next：下一个度分裂 ——
    fn deg_next(&mut self) -> bool {
        // 找最右边的度 >1
        let mut i = self.ncpt - 1;
        while i > 0 {
            if cdeg_!(self, i) > 1 {
                break;
            }
            i -= 1;
        }
        if i > 0 {
            // 找到了。若左边还有同度分量，找到最左边的并减 same
            let mut j = i - 1;
            if j > 0 && cdeg_!(self, j) == cdeg_!(self, i) {
                while j > 1 && cdeg_!(self, j - 1) == cdeg_!(self, i) {
                    j -= 1;
                }
                same!(self, j, 0) -= 1;
                for k in 1..self.ncls {
                    same!(self, j, k) = same!(self, j, 0);
                }
            }

            // 当前度 -1
            cdeg_!(self, i) -= 1;
            self.cpt_vdim[i as usize] -= 2 * (self.xdim + 2);
            self.cpt_tm[i as usize] = 0;

            // 分配剩余度
            let mut d = self.ncpt - i; // 剩余度
            let k = d / cdeg_!(self, i) + 1;
            same!(self, i, 0) = k;
            for j in 1..self.ncls {
                same!(self, i, j) = same!(self, i, 0);
            }
            i += 1;

            for _ in 2..=k {
                cdeg_!(self, i) = cdeg_!(self, i - 1);
                d -= cdeg_!(self, i);
                self.cpt_vdim[i as usize] = 2 * ((2 + self.xdim) * cdeg_!(self, i) + self.xdim - 2);
                self.cpt_tm[i as usize] = 0;
                i += 1;
            }

            // 余数
            if d != 0 {
                cdeg_!(self, i) = d;
                self.cpt_vdim[i as usize] = 2 * ((2 + self.xdim) * d + self.xdim - 2);
                self.cpt_tm[i as usize] = 0;
                same!(self, i, 0) = 1;
                for k in 1..self.ncls {
                    same!(self, i, k) = same!(self, i, 0);
                }
                i += 1;
            }

            self.ncpt = i;
            return true;
        }

        // 所有外部度都是 1：内部分量度 -1
        if cdeg_!(self, 0) == 0 {
            return false;
        }
        cdeg_!(self, 0) -= 1;
        cdeg_!(self, 1) = self.ncpt;
        self.cpt_vdim[1] = 2 * ((2 + self.xdim) * self.ncpt + self.xdim - 2);
        self.cpt_tm[1] = 0;
        self.ncpt = 2;
        for k in 0..self.ncls {
            same!(self, 1, k) = 1;
        }
        true
    }

    // —— same_next：same 块内部的下一个分裂（亏格与类共用）——
    fn same_next(&mut self, m: i32, n: i32) -> bool {
        if value!(self, m, n) == 0 {
            return false;
        }
        let big_n = same!(self, m, n - 1); // 可用分量数
        if big_n == 1 {
            return false;
        }

        let mut i = m + big_n - 1;
        let mut c = value!(self, i, n);
        let mut j = 2;
        while j <= big_n {
            i -= 1;
            let k = value!(self, i, n);
            c += k;
            if n != 1 {
                self.coeff *= factorial(value!(self, i + 1, n) as i64);
            }
            self.cpt_vdim[(i + 1) as usize] -= value!(self, i + 1, n) * self.cls_vdim[n as usize];
            if n != 1 {
                self.cpt_tm[(i + 1) as usize] -=
                    value!(self, i + 1, n) * self.cls[n as usize].mult as i32;
            }
            value!(self, i + 1, n) = 0;
            // 余量够分才能减小 k
            if k > 1 && (c + k - 2) / (k - 1) <= j {
                self.same_next_found(i, j, c, n, big_n);
                return true;
            }
            j += 1;
        }

        // 无更多分裂：复位
        if n != 1 {
            self.coeff *= factorial(value!(self, i, n) as i64) / factorial(c as i64);
        }
        self.cpt_vdim[i as usize] += (c - value!(self, i, n)) * self.cls_vdim[n as usize];
        if n != 1 {
            self.cpt_tm[i as usize] += (c - value!(self, i, n)) * self.cls[n as usize].mult as i32;
        }
        value!(self, i, n) = c;
        for k in n..self.ncls {
            same!(self, i, k) = 1;
            same!(self, i + 1, k) = big_n - 1;
        }
        false
    }

    // same_next 的 FOUND 部分（拆出以缓解借用冲突；变量语义同原版）
    fn same_next_found(&mut self, mut i: i32, mut j: i32, mut c: i32, n: i32, big_n: i32) {
        // 减小当前值前，若左边还有同值，减对应 same
        let mut kk = 0;
        while kk < big_n - j {
            if value!(self, i - kk - 1, n) != value!(self, i, n) {
                break;
            }
            kk += 1;
        }
        if kk != 0 {
            same!(self, i - kk, n) -= 1;
            for l in n + 1..self.ncls {
                same!(self, i - kk, l) = same!(self, i - kk, n);
            }
        }

        // 减小当前值
        if n != 1 {
            self.coeff *= Q::from(value!(self, i, n) as i64);
        }
        self.cpt_vdim[i as usize] -= self.cls_vdim[n as usize];
        if n != 1 {
            self.cpt_tm[i as usize] -= self.cls[n as usize].mult as i32;
        }
        value!(self, i, n) -= 1;

        // 写入对应 same 值
        let k = c / value!(self, i, n);
        same!(self, i, n) = k;
        c -= value!(self, i, n);
        for l in n + 1..self.ncls {
            same!(self, i, l) = k;
        }

        // 当前值重复写
        for _ in 1..k {
            if n != 1 {
                self.coeff /= factorial(value!(self, i, n) as i64);
            }
            self.cpt_vdim[(i + 1) as usize] += value!(self, i, n) * self.cls_vdim[n as usize];
            if n != 1 {
                self.cpt_tm[(i + 1) as usize] +=
                    value!(self, i, n) * self.cls[n as usize].mult as i32;
            }
            value!(self, i + 1, n) = value!(self, i, n);
            i += 1;
            j -= 1;
            c -= value!(self, i, n);
        }

        // 写余数
        if c != 0 {
            j -= 1;
            if n != 1 {
                self.coeff /= factorial(c as i64);
            }
            i += 1;
            self.cpt_vdim[i as usize] += c * self.cls_vdim[n as usize];
            if n != 1 {
                self.cpt_tm[i as usize] += c * self.cls[n as usize].mult as i32;
            }
            value!(self, i, n) = c;
            for l in n..self.ncls {
                same!(self, i, l) = 1;
            }
        }

        // 还剩 j 个零要写
        j -= 1;
        if j > 0 {
            i += 1;
            for l in n..self.ncls {
                same!(self, i, l) = j;
            }
            for _ in 1..=j {
                value!(self, i, n) = 0;
                i += 1;
            }
        }
    }

    // —— group_next：类在 same 组之间的下一个分配 ——
    fn group_next(&mut self, n: i32) -> bool {
        let mut m: i32 = 0;

        // 找最左非零项
        while value!(self, m, n) == 0 {
            m += same!(self, m, n - 1);
            if m >= self.ncpt {
                return false;
            }
        }

        let m2 = m + same!(self, m, n - 1);
        if m2 >= self.ncpt {
            // 循环结束：恢复初始状态
            if m == 0 {
                return false;
            }
            value!(self, 0, n) = value!(self, m, n);
            value!(self, m, n) = 0;
            self.cpt_vdim[m as usize] -= value!(self, 0, n) * self.cls_vdim[n as usize];
            self.cpt_vdim[0] += value!(self, 0, n) * self.cls_vdim[n as usize];
            if n != 1 {
                self.cpt_tm[m as usize] -= value!(self, 0, n) * self.cls[n as usize].mult as i32;
                self.cpt_tm[0] += value!(self, 0, n) * self.cls[n as usize].mult as i32;
            }
            for i in n..self.ncls {
                same!(self, m, i) = same!(self, m, n - 1);
            }
            return false;
        }

        // 右移一个类
        let was_zero = value!(self, m2, n) == 0;
        value!(self, m2, n) += 1;
        if was_zero {
            for i in n..self.ncls {
                same!(self, m2, i) = 1;
            }
            if same!(self, m2, n - 1) > 1 {
                for i in n..self.ncls {
                    same!(self, m2 + 1, i) = same!(self, m2, n - 1) - 1;
                }
            }
        }
        if n != 1 {
            self.coeff *= Q::from(value!(self, m, n) as i64);
            self.coeff /= Q::from(value!(self, m2, n) as i64);
            self.cpt_tm[m as usize] -= self.cls[n as usize].mult as i32;
            self.cpt_tm[m2 as usize] += self.cls[n as usize].mult as i32;
        }
        value!(self, m, n) -= 1;
        self.cpt_vdim[m as usize] -= self.cls_vdim[n as usize];
        self.cpt_vdim[m2 as usize] += self.cls_vdim[n as usize];

        // 余类移到最左
        if m != 0 && value!(self, m, n) != 0 {
            value!(self, 0, n) = value!(self, m, n);
            value!(self, m, n) = 0;
            self.cpt_vdim[m as usize] -= value!(self, 0, n) * self.cls_vdim[n as usize];
            if n != 1 {
                self.cpt_tm[m as usize] -= value!(self, 0, n) * self.cls[n as usize].mult as i32;
            }
        }

        // 更新当前位置 same
        if value!(self, m, n) == 0 {
            for i in n..self.ncls {
                same!(self, m, i) = same!(self, m, n - 1);
            }
        }

        true
    }

    // —— loop_next：下一个环路节点配置 ——
    fn loop_next(&mut self) -> bool {
        let mut j = 1i32; // 当前 same 块起点
        let mut k = 1 + same!(self, 1, 1); // 下一块起点
        let mut i = 1i32;
        while i < self.ncpt {
            if i >= k {
                j = k;
                k += same!(self, i, 1);
            }
            cloop_!(self, i) += 1;
            if cloop_!(self, i) <= cgenus_!(self, i) {
                // FOUND
                self.cpt_vdim[i as usize] -= self.cls_vdim[1];
                for k2 in j..i {
                    self.cpt_vdim[k2 as usize] -=
                        (cloop_!(self, k2) - cloop_!(self, i)) * self.cls_vdim[1];
                    cloop_!(self, k2) = cloop_!(self, i);
                }
                if same!(self, i, 2) > 1 {
                    for l in 2..self.ncls {
                        same!(self, i + 1, l) = same!(self, i, 2) - 1;
                    }
                }
                for l in 2..self.ncls {
                    same!(self, j, l) = i - j + 1;
                }
                return true;
            }
            // 溢出
            self.cpt_vdim[i as usize] += (cloop_!(self, i) - 1) * self.cls_vdim[1];
            cloop_!(self, i) = 0;
            i += 1;
        }
        false
    }

    // —— nodemult_init：检查分量 ci 的节点配置合法性并初始化重数 ——
    fn nodemult_init(&mut self, ci: i32) -> bool {
        // nodeclass_init 开头的特例（见下）
        if self.noclass && cloop_!(self, ci) == 0 {
            return true;
        }

        // 原始类必须成对
        if self.xdeg > 1 {
            // 先查非节点标记点上的原始类
            let mut j = 0i32;
            if self.cpt_prim[0] != 0 {
                for i in 3..self.ncls {
                    let mut k;
                    if (self.cls[i as usize].h as i32) < j {
                        j = self.cls[i as usize].h as i32;
                        k = 0;
                        for l in i..self.ncls {
                            if (self.cls[l as usize].h as i32) == j {
                                k += value!(self, ci, l);
                            }
                        }
                        for l in 0..=cloop_!(self, ci) {
                            if node_class!(self, ci, l) == j {
                                k += 1;
                            }
                        }
                        if k & 1 != 0 {
                            return false;
                        }
                    }
                }
            }

            // 只出现在节点上的原始类（注意：原版此处用 cloop(i)，直译保留）
            for i in 0..=cloop_!(self, ci) {
                if node_class!(self, ci, i) < j {
                    j = node_class!(self, ci, i);
                    let mut k = 1i32;
                    let mut l = i + 1;
                    while l <= cloop_!(self, i) {
                        if node_class!(self, ci, l) == j {
                            k += 1;
                        }
                        l += 1;
                    }
                    if k & 1 != 0 {
                        return false;
                    }
                }
            }
        }

        // 节点总重数
        let mut vd = self.cpt_vdim[ci as usize];
        for j in 0..=cloop_!(self, ci) {
            vd -= if node_class!(self, ci, j) < 0 {
                self.xdim
            } else {
                2 * node_class!(self, ci, j)
            };
        }
        // 原始类已保证成对，虚维数必为偶
        vd = vd / 2 + cloop_!(self, ci) + 1;
        if self.cpt_tm[ci as usize] + vd > self.xdeg * cdeg_!(self, ci) {
            return false; // 总重数超限
        }

        // 首重数最大，其余为 1
        node_mult!(self, ci, 0) = vd - cloop_!(self, ci);
        if node_mult!(self, ci, 0) < 1 {
            return false;
        }
        for j in 1..=cloop_!(self, ci) {
            node_mult!(self, ci, j) = 1;
        }

        true
    }

    // —— nodemult_next：分量 ci 的下一个重数分裂 ——
    fn nodemult_next(&mut self, ci: i32) -> bool {
        if cloop_!(self, ci) == 0 {
            return false; // 单节点无自由度
        }

        // 先在 node_same 块内部变更
        let mut ni = 0i32;
        while ni <= cloop_!(self, ci) {
            let big_n = node_same!(self, ci, ni);
            if big_n == 1 {
                ni += big_n;
                continue;
            }

            let mut c = node_mult!(self, ci, ni + big_n - 1) - 1;
            let mut i = ni + big_n - 2;
            while i >= ni {
                let k = node_mult!(self, ci, i);
                c += k - 1;
                node_mult!(self, ci, i + 1) = 1;
                if k > 2 && (c + (k - 3)) / (k - 2) <= ni + big_n - i {
                    // FOUND：减小重数并分配
                    node_mult!(self, ci, i) -= 1;
                    let k = k - 1; // 原：c -= (--k)-1
                    c -= k - 1;
                    while c >= k - 1 {
                        i += 1;
                        node_mult!(self, ci, i) = k;
                        c -= k - 1;
                    }
                    if c != 0 {
                        i += 1;
                        node_mult!(self, ci, i) = c + 1;
                    }
                    return true;
                }
                i -= 1;
            }

            // 本块无更多分裂：复位，试下一块
            node_mult!(self, ci, ni) = c + 1;
            ni += big_n;
        }

        // 块内分布已尽：改变块间分配。找最左非 1 项
        let mut ni = 0i32;
        while node_mult!(self, ci, ni) == 1 {
            ni += node_same!(self, ci, ni);
            if ni > cloop_!(self, ci) {
                return false;
            }
        }

        let i = ni + node_same!(self, ci, ni);
        if i > cloop_!(self, ci) {
            return false;
        }

        node_mult!(self, ci, i) += 1;
        node_mult!(self, ci, ni) -= 1;

        if ni == 0 {
            return true;
        }
        node_mult!(self, ci, 0) += node_mult!(self, ci, ni) - 1;
        node_mult!(self, ci, ni) = 1;
        true
    }

    // —— nodeclass_init：分量 ci 的节点类初始化 ——
    fn nodeclass_init(&mut self, ci: i32) -> bool {
        // 特例：无外部类且单节点（同时设好重数，nodemult_init 直接通过）
        if self.noclass && cloop_!(self, ci) == 0 {
            let mut j = self.xdeg * cdeg_!(self, ci);
            let mut i = self.cpt_vdim[ci as usize] / 2 - (j - 1);
            if i > self.xdim {
                return false;
            }
            if i < 0 {
                j += i;
                i = 0;
                if j < 1 {
                    return false;
                }
            }
            node_class!(self, ci, 0) = i;
            node_mult!(self, ci, 0) = j;
            node_same!(self, ci, 0) = 1;
            return true;
        }

        // 所有节点取最小类
        for j in 0..=cloop_!(self, ci) {
            node_class!(self, ci, j) = self.node_minclass;
        }
        node_same!(self, ci, 0) = cloop_!(self, ci) + 1;
        true
    }

    // —— nodeclass_next：分量 ci 的下一个节点类配置 ——
    fn nodeclass_next(&mut self, ci: i32) -> bool {
        // 特例：无外部类且单节点
        if self.noclass && cloop_!(self, ci) == 0 {
            if node_class!(self, ci, 0) == self.xdim || node_mult!(self, ci, 0) == 1 {
                return false;
            }
            node_class!(self, ci, 0) += 1;
            node_mult!(self, ci, 0) -= 1;
            return true;
        }

        for ni in 0..=cloop_!(self, ci) {
            let mut raised = false;

            // 还有更大 H 幂：乘 H
            if node_class!(self, ci, ni) >= 0 && node_class!(self, ci, ni) < self.xdim {
                node_class!(self, ci, ni) += 1;
                raised = true;
            } else if node_class!(self, ci, ni) == self.xdim
                && self.xdeg > 1
                && self.xdim >= 2 * self.node_minclass
                && (cloop_!(self, ci) != 0 || self.cpt_prim[0] != 0)
            {
                // 切到原始类
                node_class!(self, ci, ni) = -1;
                raised = true;
            } else if node_class!(self, ci, ni) < 0
                && (node_class!(self, ci, ni) >= self.cpt_prim[ci as usize]
                    || node_same!(self, ci, ni) > 1)
            {
                // 下一个原始类
                node_class!(self, ci, ni) -= 1;
                raised = true;
            }

            if raised {
                // 更新后续节点的 node_same 与前驱节点的类
                if node_same!(self, ci, ni) > 1 {
                    node_same!(self, ci, ni + 1) = node_same!(self, ci, ni) - 1;
                }
                for i in 0..ni {
                    node_class!(self, ci, i) = node_class!(self, ci, ni);
                }
                node_same!(self, ci, 0) = ni + 1;
                return true;
            }
        }
        false
    }

    // —— node_init：全部分量的节点初始化 ——
    fn node_init(&mut self) -> bool {
        let mut i = 1i32;
        while i < self.ncpt {
            if !self.nodeclass_init(i) {
                return false;
            }
            while !self.nodemult_init(i) {
                if !self.nodeclass_next(i) {
                    return false;
                }
            }

            // 拷贝到 same 块
            let block = same!(self, i, self.ncls - 1);
            for j in i + 1..i + block {
                for k in 0..=cloop_!(self, j) {
                    node_class!(self, j, k) = node_class!(self, i, k);
                    node_mult!(self, j, k) = node_mult!(self, i, k);
                    node_same!(self, j, k) = node_same!(self, i, k);
                }
            }

            // 下一 same 块（写 ncls 与 ncls+1 两层）
            let s = same!(self, i, self.ncls - 1);
            same!(self, i, self.ncls + 1) = s;
            same!(self, i, self.ncls) = s;
            i += s;
        }
        true
    }

    // —— node_next：下一个节点（类+重数）配置 ——
    fn node_next(&mut self) -> bool {
        // 先试重数配置
        let mut ci = 1i32;
        let mut cj = 1i32;
        while ci < self.ncpt {
            if self.nodemult_next(ci) {
                for i in cj..ci {
                    for j in 0..=cloop_!(self, i) {
                        node_mult!(self, i, j) = node_mult!(self, ci, j);
                    }
                }
                if same!(self, ci, self.ncls + 1) > 1 {
                    same!(self, ci + 1, self.ncls + 1) = same!(self, ci, self.ncls + 1) - 1;
                }
                same!(self, cj, self.ncls + 1) = ci - cj + 1;
                return true;
            }

            ci += 1;
            if ci == cj + same!(self, cj, self.ncls) {
                // 越过一个 same 块：复位
                self.nodemult_init(cj);
                for i in cj + 1..ci {
                    for j in 0..=cloop_!(self, i) {
                        node_mult!(self, i, j) = node_mult!(self, cj, j);
                    }
                }
                same!(self, cj, self.ncls + 1) = same!(self, cj, self.ncls);
                cj = ci;
            }
        }

        // 重数配置已尽：试类配置
        let mut ci = 1i32;
        let mut cj = 1i32;
        while ci < self.ncpt {
            loop {
                if !self.nodeclass_next(ci) {
                    break;
                }
                if self.nodemult_init(ci) {
                    for i in cj..ci {
                        for j in 0..=cloop_!(self, i) {
                            node_mult!(self, i, j) = node_mult!(self, ci, j);
                            node_class!(self, i, j) = node_class!(self, ci, j);
                            node_same!(self, i, j) = node_same!(self, ci, j);
                        }
                    }
                    if same!(self, ci, self.ncls) > 1 {
                        let v = same!(self, ci, self.ncls) - 1;
                        same!(self, ci + 1, self.ncls + 1) = v;
                        same!(self, ci + 1, self.ncls) = v;
                    }
                    same!(self, cj, self.ncls + 1) = ci - cj + 1;
                    same!(self, cj, self.ncls) = ci - cj + 1;
                    // 更新 cpt_prim
                    for i in (1..ci).rev() {
                        self.cpt_prim[i as usize] = self.cpt_prim[i as usize + 1];
                        if node_class!(self, i + 1, 0) < self.cpt_prim[i as usize] {
                            self.cpt_prim[i as usize] = node_class!(self, i + 1, 0);
                        }
                    }
                    return true;
                }
            }

            ci += 1;
            if ci == cj + same!(self, cj, self.ncls - 1) {
                if !self.nodeclass_init(cj) {
                    return false;
                }
                while !self.nodemult_init(cj) {
                    if !self.nodeclass_next(cj) {
                        return false;
                    }
                }
                for i in cj + 1..ci {
                    for j in 0..=cloop_!(self, i) {
                        node_mult!(self, i, j) = node_mult!(self, cj, j);
                        node_class!(self, i, j) = node_class!(self, cj, j);
                        node_same!(self, i, j) = node_same!(self, cj, j);
                    }
                }
                let s = same!(self, cj, self.ncls - 1);
                same!(self, cj, self.ncls + 1) = s;
                same!(self, cj, self.ncls) = s;
                cj = ci;
            }
        }
        false
    }

    // —— 主迭代：下一个相对分裂 ——
    // 原版 relinv::next 的直译：推进状态机，返回是否还有下一个分裂（不是迭代器）
    #[allow(clippy::should_implement_trait)]
    pub fn next(&mut self) -> bool {
        if self.node_minclass > self.xdim {
            return false;
        }
        if crate::env_flag("SC_GROWI_DEBUG") {
            use std::sync::atomic::{AtomicU64, Ordering};
            static N: AtomicU64 = AtomicU64::new(0);
            let n = N.fetch_add(1, Ordering::Relaxed);
            if n % 1_000_000 == 0 {
                eprintln!(
                    "[relinv #{n}] ncpt={} cdeg0={} ncls={} coeff={}",
                    self.ncpt,
                    cdeg_!(self, 0),
                    self.ncls,
                    self.coeff
                );
            }
        }

        use St::*;
        let mut st = Begin;
        loop {
            st = match st {
                Begin => {
                    if self.ncpt == 1 {
                        Degree
                    } else {
                        NodeNextTry
                    }
                }
                NodeNextTry => {
                    if self.node_next() {
                        DoneCheck
                    } else {
                        Again
                    }
                }
                Again => {
                    if self.advance_classes() {
                        Ininod
                    } else {
                        Degree
                    }
                }
                Degree => {
                    if self.deg_next() {
                        Ininod
                    } else {
                        return false;
                    }
                }
                Ininod => {
                    self.fix_big_classes();
                    if self.node_init() { DoneCheck } else { Again }
                }
                DoneCheck => {
                    if self.mult_cond() {
                        return true;
                    }
                    Begin
                }
            };
        }
    }

    fn advance_classes(&mut self) -> bool {
        // 类分裂（n 从高到低）
        for n in (3..self.ncls).rev() {
            let mut m = 1i32;
            while m < self.ncpt {
                if self.same_next(m, n) {
                    return true;
                }
                m += same!(self, m, n - 1);
            }
            if self.group_next(n) {
                return true;
            }
        }
        // 亏格/环路分裂（总亏格非零时）
        if self.totalgenus != 0 {
            if self.loop_next() {
                return true;
            }
            let mut m = 0i32;
            while m < self.ncpt {
                if self.same_next(m, 1) {
                    return true;
                }
                m += same!(self, m, 0);
            }
            if self.group_next(1) {
                return true;
            }
        }
        false
    }

    /// ININOD：把内部分量上 H 幂过大的类挪到分量 1。
    fn fix_big_classes(&mut self) {
        let mut changed = false;
        for i in 3..self.ncls {
            if (self.cls[i as usize].h as i32) > self.xdim && value!(self, 0, i) != 0 {
                self.coeff *= factorial(value!(self, 0, i) as i64);
                self.cpt_vdim[1] += value!(self, 0, i) * self.cls_vdim[i as usize];
                self.cpt_tm[1] += value!(self, 0, i) * self.cls[i as usize].mult as i32;
                value!(self, 1, i) += value!(self, 0, i);
                value!(self, 0, i) = 0;
                self.coeff /= factorial(value!(self, 1, i) as i64);
                changed = true;
            }
            if changed && same!(self, 1, i) > 1 {
                same!(self, 2, i) = same!(self, 1, i) - 1;
                same!(self, 1, i) = 1;
            }
        }
    }

    /// 取分量 n 的不变量。n=0 为内部分量：节点带负重数、内部点带正
    /// 重数（更高亏格相对分裂需要），调用方 compute 前必须 delmult()。
    pub fn cpt(&self, n: i32) -> Inv {
        let mut pts: Vec<Pdata> = Vec::new();

        if n != 0 {
            // 外部分量：先节点，后其它类
            for i in 0..=cloop_!(self, n) {
                pts.push(Pdata::new(
                    node_class!(self, n, i) as i16,
                    0,
                    node_mult!(self, n, i) as i16,
                ));
            }
            for i in 3..self.ncls {
                for _ in 0..value!(self, n, i) {
                    pts.push(self.cls[i as usize]);
                }
            }
        } else {
            // 内部分量：先节点（对角类，H 幂取补、重数取负）
            for i in 1..self.ncpt {
                for j in 0..=cloop_!(self, i) {
                    let hp = if node_class!(self, i, j) >= 0 {
                        self.xdim - node_class!(self, i, j)
                    } else {
                        node_class!(self, i, j)
                    };
                    pts.push(Pdata::new(hp as i16, 0, -(node_mult!(self, i, j) as i16)));
                }
            }
            // 跳过类
            if self.skip {
                pts.push(self.skipclass);
            }
            // 其它类
            for i in 3..self.ncls {
                for _ in 0..value!(self, 0, i) {
                    pts.push(self.cls[i as usize]);
                }
            }
        }

        Inv::new(
            self.xdim as i16,
            self.xdeg as i16,
            (cgenus_!(self, n) - cloop_!(self, n)) as i16,
            cdeg_!(self, n) as i16,
            pts,
        )
    }

    /// 当前分裂的数值系数（组合因子 + 对角分裂 + 重数因子 + 原始类计数）。
    pub fn cof(&self) -> Q {
        let mut result = self.coeff.clone();

        // 相同分量因子
        let mut i = 1i32;
        while i < self.ncpt {
            if same!(self, i, self.ncls + 1) > 1 {
                result /= factorial(same!(self, i, self.ncls + 1) as i64);
            }
            i += same!(self, i, self.ncls + 1);
        }

        // 相同节点因子
        if self.totalgenus != 0 {
            for i in 1..self.ncpt {
                let mut j = 1i32;
                let mut c = 0i32;
                for k in 0..=cloop_!(self, i) {
                    j -= 1;
                    if j == 0 {
                        j = node_same!(self, i, k);
                        c = 0;
                    }
                    if k > 0 && node_mult!(self, i, k - 1) != node_mult!(self, i, k) {
                        c = 0;
                    }
                    c += 1;
                    if c > 1 {
                        result /= Q::from(c as i64);
                    }
                }
            }
        }

        // 对角分裂与重数因子
        for i in 1..self.ncpt {
            for j in 0..=cloop_!(self, i) {
                if node_class!(self, i, j) >= 0 {
                    result /= Q::from(self.xdeg as i64);
                }
                if node_mult!(self, i, j) > 1 {
                    result *= Q::from(node_mult!(self, i, j) as i64);
                }
            }
        }

        // 原始类个数因子
        if self.xdeg > 1 && (self.totalgenus != 0 || self.cpt_prim[0] != 0) {
            let mut j = self.cpt_prim[1];
            if node_class!(self, 1, 0) < j {
                j = node_class!(self, 1, 0);
            }
            if self.prims.clone() + Q::from(j as i64) < Q::ZERO {
                return Q::ZERO; // 用了超过存在的原始类个数
            }
            for i in self.cpt_prim[0]..j {
                result *= self.prims.clone() + Q::from(i as i64);
            }
            for i in 1..self.ncpt {
                result /= factorial((self.cpt_prim[i as usize] - j) as i64);
                j = self.cpt_prim[i as usize];
            }
        }

        result
    }

    /// 相对 GW 理论的"重数条件"：外部重数之和 = 内部重数之和 + 交点数？
    /// 同时检查内部分量无超限 H 幂。
    pub fn mult_cond(&self) -> bool {
        let mut n = self.xdeg * cdeg_!(self, 0);

        // 内部分量 H 幂
        for i in 3..self.ncls {
            if (self.cls[i as usize].h as i32) > self.xdim && value!(self, 0, i) != 0 {
                return false;
            }
        }

        // 外部重数
        for i in 1..self.ncpt {
            for j in 0..=cloop_!(self, i) {
                n += node_mult!(self, i, j);
            }
        }

        // 内部重数
        for i in 3..self.ncls {
            if self.cls[i as usize].mult != 0 {
                n -= value!(self, 0, i) * self.cls[i as usize].mult as i32;
            }
        }
        if self.skip {
            n -= self.skipclass.mult as i32;
        }

        n == 0
    }
}

/// next() 的控制流状态（直译原版的 BEGIN/AGAIN/DEGREE/ININOD/DONE 标签）。
enum St {
    Begin,
    NodeNextTry,
    Again,
    Degree,
    Ininod,
    DoneCheck,
}

/// 整数底有理数幂（指数可为负）。
fn q_pow(x: Q, n: i64) -> Q {
    let neg = n < 0;
    let mut e = n.unsigned_abs();
    let mut base = x;
    let mut acc = Q::ONE;
    while e > 0 {
        if e & 1 == 1 {
            acc = &acc * &base;
        }
        base = &base * &base;
        e >>= 1;
    }
    if neg { Q::ONE / &acc } else { acc }
}
