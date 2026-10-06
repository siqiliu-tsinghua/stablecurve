//! 不变量 hash-consing 仓——取代原版"排序 + 位打包 + gdbm 键"三件套
//! （REDESIGN §2.2）。
//!
//! 规范形 [`Inv`] 全局唯一化为 [`InvId`]（u32 句柄）：相等性退化为整数比较，
//! 克隆免费，记忆化缓存（M1 的 L1/L2）直接以 id 为键。这是整个重做中单项
//! 收益最大的结构改动。

use std::collections::HashMap;

use crate::inv::Inv;

/// 不变量的进程内唯一句柄。仅在 [`InvStore`] 内有意义。
#[derive(Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Debug)]
pub struct InvId(u32);

impl InvId {
    #[inline]
    pub fn idx(self) -> usize {
        self.0 as usize
    }
}

/// hash-consing 仓。`intern` 会先把不变量排成规范形再登记，因此仓内
/// 永远只有规范形——两个语义相同的不变量必得同一 id。
#[derive(Default)]
pub struct InvStore {
    map: HashMap<Inv, InvId>,
    slots: Vec<Inv>,
}

impl InvStore {
    pub fn new() -> Self {
        InvStore::default()
    }

    /// 登记一个不变量（自动规范形化），返回其唯一 id。
    pub fn intern(&mut self, inv: Inv) -> InvId {
        if let Some(&id) = self.map.get(&inv) {
            return id;
        }
        let mut inv = inv;
        inv.canonicalize();
        if let Some(&id) = self.map.get(&inv) {
            return id;
        }
        let id = InvId(self.slots.len() as u32);
        self.slots.push(inv.clone());
        self.map.insert(inv, id);
        id
    }

    /// 按 id 取回规范形不变量。
    ///
    /// # Panics
    /// id 非本仓签发时 panic（内部一致性错误）。
    pub fn get(&self, id: InvId) -> &Inv {
        &self.slots[id.idx()]
    }

    pub fn len(&self) -> usize {
        self.slots.len()
    }

    pub fn is_empty(&self) -> bool {
        self.slots.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::inv::Pdata;

    #[test]
    fn same_invariant_same_id() {
        let mut store = InvStore::new();
        let a = Inv::new(2, 1, 0, 1, [Pdata::new(2, 0, 0); 2]);
        let b = a.clone();
        assert_eq!(store.intern(a), store.intern(b));
        assert_eq!(store.len(), 1);
    }

    #[test]
    fn canonical_forms_converge() {
        // 点顺序不同、语义相同 → 同一 id
        let mut store = InvStore::new();
        let a = Inv::new(3, 1, 0, 2, [Pdata::new(1, 0, 0), Pdata::new(2, 0, 0)]);
        let b = Inv::new(3, 1, 0, 2, [Pdata::new(2, 0, 0), Pdata::new(1, 0, 0)]);
        let ia = store.intern(a);
        let ib = store.intern(b);
        assert_eq!(ia, ib);
        // 仓里存的是规范形（H1 在前）
        assert_eq!(store.get(ia).pts[0].h, 1);
    }

    #[test]
    fn different_invariants_different_ids() {
        let mut store = InvStore::new();
        let a = Inv::new(2, 1, 0, 1, [Pdata::new(2, 0, 0); 2]);
        let b = Inv::new(2, 1, 0, 2, [Pdata::new(2, 0, 0); 2]);
        assert_ne!(store.intern(a), store.intern(b));
        assert_eq!(store.len(), 2);
    }
}
