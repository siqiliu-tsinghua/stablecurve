//! Backend-only reads of the documented pre-source SQLite representations.
//! Old tables remain intact. Only legacy KLP roots require a lazy, bounded
//! import because their historical keys preceded input normalization.

use super::IntegralKey;
use super::codec::{encode_key_into, encode_unranked_key_into};
use rusqlite::{Connection, OptionalExtension, params};
use sc_num::{Q, decode_rational, encode_rational_into};

const LEGACY: &[u8] = b"sc-faber 0.1.0|";
const DESCENDANT: &[u8] = b"psi-grr|";
const OLD_DESCENDANT: &[u8] = b"psi-normalized-grr|";
const KLP: &[u8] = b"sc-faber 0.1.0|klp|";
const KLP_END: &[u8] = b"sc-faber 0.1.0|klp}";
const IMPORT_ROWS: usize = 512;

pub(super) struct LegacyReader {
    has_values: bool,
    klp_imported: bool,
}

fn sql_error(error: rusqlite::Error) -> String {
    format!("legacy FABER cache: {error}")
}

fn put(mut value: u64, bytes: &mut Vec<u8>) {
    while value >= 128 {
        bytes.push((value as u8 & 127) | 128);
        value >>= 7;
    }
    bytes.push(value as u8);
}

fn length(value: usize) -> Result<u32, String> {
    value
        .try_into()
        .map_err(|_| "legacy key length exceeds u32".into())
}

fn old_key(kind: &str, nums: &[u32]) -> Vec<u8> {
    let mut bytes = LEGACY.to_vec();
    bytes.extend_from_slice(kind.as_bytes());
    bytes.push(b'|');
    for &number in nums {
        put(u64::from(number), &mut bytes);
    }
    bytes
}

fn read_numbers(mut bytes: &[u8]) -> Result<Vec<u32>, String> {
    let mut numbers = Vec::new();
    while !bytes.is_empty() {
        let mut value = 0u32;
        let mut complete = false;
        for index in 0..5 {
            let (&byte, rest) = bytes.split_first().ok_or("truncated legacy key")?;
            bytes = rest;
            if index == 4 && byte > 15 {
                return Err("legacy key integer exceeds u32".into());
            }
            value |= u32::from(byte & 127) << (7 * index);
            if byte < 128 {
                complete = true;
                break;
            }
        }
        if !complete {
            return Err("legacy key integer exceeds u32".into());
        }
        numbers.push(value);
    }
    Ok(numbers)
}

fn legacy_identity(bytes: &[u8]) -> Result<Option<IntegralKey>, String> {
    let Some(body) = bytes.strip_prefix(LEGACY) else {
        return Ok(None);
    };
    let boundary = body
        .iter()
        .position(|&byte| byte == b'|')
        .ok_or("legacy key has no kind boundary")?;
    let kind = std::str::from_utf8(&body[..boundary]).map_err(|_| "legacy kind is not UTF-8")?;
    Ok(Some(IntegralKey::Legacy {
        kind: kind.to_owned(),
        nums: read_numbers(&body[boundary + 1..])?,
    }))
}

/// Convert the old sign/u32-LE-numerator-length/rest-denominator envelope to
/// wire version 1, then use the shared exact decoder. No decimal intermediate.
fn old_rational(bytes: &[u8]) -> Result<Q, String> {
    if bytes.len() < 6 || bytes[0] > 1 {
        return Err("invalid legacy rational envelope".into());
    }
    let numerator = u32::from_le_bytes(bytes[1..5].try_into().unwrap()) as usize;
    let end = 5usize
        .checked_add(numerator)
        .filter(|&end| end < bytes.len())
        .ok_or("invalid legacy rational numerator length")?;
    let mut wire = vec![1, bytes[0]];
    put(numerator as u64, &mut wire);
    wire.extend_from_slice(&bytes[5..end]);
    put((bytes.len() - end) as u64, &mut wire);
    wire.extend_from_slice(&bytes[end..]);
    decode_rational(&wire)
}

fn current_value(connection: &Connection, key: &IntegralKey) -> Result<Option<Q>, String> {
    let mut encoded = Vec::new();
    encode_key_into(key, &mut encoded);
    if let Some(value) = typed_value(connection, &encoded)? {
        return Ok(Some(value));
    }
    unranked_value(connection, key)
}

fn current_value_unless_tried(
    connection: &Connection,
    key: &IntegralKey,
    typed_already_tried: Option<&IntegralKey>,
) -> Result<Option<Q>, String> {
    if typed_already_tried == Some(key) {
        return Ok(None);
    }
    current_value(connection, key)
}

fn typed_value(connection: &Connection, encoded: &[u8]) -> Result<Option<Q>, String> {
    let value: Option<Vec<u8>> = connection
        .prepare_cached("SELECT v FROM faber_integrals WHERE k=?1")
        .map_err(sql_error)?
        .query_row([encoded], |row| row.get(0))
        .optional()
        .map_err(sql_error)?;
    value.map(|bytes| decode_rational(&bytes)).transpose()
}

/// The old typed descendant retained the same C normalization and full
/// mathematical identity. Read its exact primary key only after a new miss.
fn unranked_value(connection: &Connection, key: &IntegralKey) -> Result<Option<Q>, String> {
    let IntegralKey::Descendant { g, d, .. } = key else {
        return Ok(None);
    };
    if 2 * i128::from(*g) - 2 + d.len() as i128 <= 0 {
        // Auxiliary unstable descendants already retain tag 1 in the current
        // encoder, so the preceding ordinary lookup checked this same key.
        return Ok(None);
    }
    let mut encoded = Vec::new();
    encode_unranked_key_into(key, &mut encoded);
    typed_value(connection, &encoded)
}

/// Exact aliases into the current typed table only. Do not call compatibility
/// entry points here: typed-to-legacy fallback already uses this reader.
fn typed_alias(
    connection: &Connection,
    kind: &str,
    nums: &[u32],
    typed_already_tried: Option<&IntegralKey>,
) -> Result<Option<Q>, String> {
    match kind {
        "st" => {
            let n = nums.len() as i128;
            let degree: i128 = nums.iter().map(|&power| i128::from(power)).sum();
            let triple_g = degree - n + 3;
            if triple_g < 0 || triple_g % 3 != 0 {
                return Ok(None);
            }
            let Ok(g) = u32::try_from(triple_g / 3) else {
                return Ok(None);
            };
            if 2 * i128::from(g) - 2 + n <= 0 {
                // In particular ST([])=1 is an auxiliary convention, not a
                // stable pure-psi integral on the inferred (1,0) space.
                return Ok(None);
            }
            let mut d = nums.to_vec();
            d.sort_unstable();
            current_value_unless_tried(
                connection,
                &IntegralKey::Descendant { g, d, ch: vec![] },
                typed_already_tried,
            )
        }
        "ic" => {
            let Some(root) = ic_identity(nums)? else {
                return Ok(None);
            };
            if let Some(value) = current_value_unless_tried(connection, &root, typed_already_tried)?
            {
                return Ok(Some(value));
            }
            let IntegralKey::Integral {
                g,
                n,
                kappa,
                psi,
                ch,
                ..
            } = &root
            else {
                unreachable!()
            };
            if kappa.is_empty() && ch.is_empty() {
                if let Some(d) = pure_root(*g, *n, psi)? {
                    return current_value_unless_tried(
                        connection,
                        &IntegralKey::Descendant {
                            g: *g,
                            d,
                            ch: vec![],
                        },
                        typed_already_tried,
                    );
                }
            }
            // An IC ch insertion is the actual Chern character. A nonempty
            // Descendant ch field contains normalized C and requires scaling,
            // which is deliberately outside this exact-alias adapter.
            Ok(None)
        }
        _ => Ok(None),
    }
}

fn ic_identity(nums: &[u32]) -> Result<Option<IntegralKey>, String> {
    if nums.len() < 5 || nums[2..5].iter().any(|count| count % 2 != 0) {
        return Err("invalid legacy IC component lengths".into());
    }
    let count = nums[2..5]
        .iter()
        .try_fold(5usize, |total, &count| total.checked_add(count as usize))
        .ok_or("legacy IC component lengths overflow")?;
    if count != nums.len() {
        return Err("invalid legacy IC component lengths".into());
    }
    let (g, n) = (nums[0], nums[1]);
    let mut offset = 5;
    let mut component = |count: u32| {
        let end = offset + count as usize;
        let pairs = nums[offset..end]
            .chunks_exact(2)
            .map(|pair| (pair[0], pair[1]))
            .collect();
        offset = end;
        pairs
    };
    let kappa = normalized(component(nums[2]))?;
    let raw_psi: Vec<(u32, u32)> = component(nums[3]);
    let ch = normalized(component(nums[4]))?;
    if raw_psi.iter().any(|&(leg, _)| leg == 0 || leg > n)
        || kappa.iter().any(|&(index, _)| index == 0)
        || ch.iter().any(|&(index, _)| index == 0 || index % 2 == 0)
        || 2 * i64::from(g) - 2 + i64::from(n) <= 0
    {
        return Ok(None);
    }
    let psi = normalized(raw_psi)?;
    let degree = kappa
        .iter()
        .chain(&ch)
        .try_fold(0u128, |degree, &(index, power)| {
            degree.checked_add(u128::from(index) * u128::from(power))
        })
        .and_then(|degree| {
            psi.iter().try_fold(degree, |degree, &(_, power)| {
                degree.checked_add(u128::from(power))
            })
        })
        .ok_or("legacy IC degree overflow")?;
    let top = 3 * i128::from(g) - 3 + i128::from(n);
    if degree != top as u128 {
        return Ok(None);
    }
    Ok(Some(IntegralKey::Integral {
        g,
        n,
        kappa,
        psi,
        lambda: vec![],
        ch,
    }))
}

impl LegacyReader {
    /// Positive-C descendants have no ST/IC/Integral aliases. Only stable old
    /// tag-1 identities and the two historical normalized-GRR namespaces can
    /// bypass the current ranked key. Current unstable tag-1 auxiliaries do
    /// not disable the hint. All probes seek existing primary-key ranges.
    pub(super) fn has_positive_descendant_aliases(
        &self,
        connection: &Connection,
    ) -> Result<bool, String> {
        let old_typed: bool = connection
            .query_row(
                "SELECT EXISTS(SELECT 1 FROM faber_integrals WHERE k>=x'010102' AND k<x'0102')
             OR EXISTS(SELECT 1 FROM faber_integrals WHERE k>=x'01010101' AND k<x'010102')
             OR EXISTS(SELECT 1 FROM faber_integrals WHERE k>=x'01010003' AND k<x'010101')",
                [],
                |row| row.get(0),
            )
            .map_err(sql_error)?;
        if old_typed || !self.has_values {
            return Ok(old_typed);
        }
        connection
            .query_row(
                "SELECT EXISTS(SELECT 1 FROM faber_values WHERE k>=?1 AND k<?2)
             OR EXISTS(SELECT 1 FROM faber_values WHERE k>=?3 AND k<?4)",
                params![
                    DESCENDANT,
                    b"psi-grr}".as_slice(),
                    OLD_DESCENDANT,
                    b"psi-normalized-grr}".as_slice()
                ],
                |row| row.get(0),
            )
            .map_err(sql_error)
    }

    pub(super) fn new(connection: &Connection) -> Result<Self, String> {
        let has_values = connection.query_row(
            "SELECT EXISTS(SELECT 1 FROM sqlite_master WHERE type='table' AND name='faber_values')",
            [], |row| row.get(0),
        ).map_err(sql_error)?;
        Ok(Self {
            has_values,
            klp_imported: false,
        })
    }

    fn blob(&self, connection: &Connection, key: &[u8]) -> Result<Option<Q>, String> {
        if !self.has_values {
            return Ok(None);
        }
        let bytes: Option<Vec<u8>> = connection
            .prepare_cached("SELECT v FROM faber_values WHERE k=?1")
            .map_err(sql_error)?
            .query_row([key], |row| row.get(0))
            .optional()
            .map_err(sql_error)?;
        bytes.map(|bytes| old_rational(&bytes)).transpose()
    }

    fn old_value(&self, connection: &Connection, key: &[u8]) -> Result<Option<Q>, String> {
        if let Some(value) = self.blob(connection, key)? {
            return Ok(Some(value));
        }
        let text: Option<String> = connection
            .prepare_cached("SELECT v FROM st WHERE k=?1")
            .map_err(sql_error)?
            .query_row([key], |row| row.get(0))
            .optional()
            .map_err(sql_error)?;
        text.map(|text| {
            text.parse::<Q>()
                .map_err(|_| "invalid legacy rational text".into())
        })
        .transpose()
    }

    #[cfg(test)]
    pub(super) fn legacy_get(
        &self,
        connection: &Connection,
        key: &[u8],
    ) -> Result<Option<Q>, String> {
        self.legacy_get_unless_tried(connection, key, None)
    }

    fn legacy_get_unless_tried(
        &self,
        connection: &Connection,
        key: &[u8],
        typed_already_tried: Option<&IntegralKey>,
    ) -> Result<Option<Q>, String> {
        if let Some(identity) = legacy_identity(key)? {
            if let Some(value) =
                current_value_unless_tried(connection, &identity, typed_already_tried)?
            {
                return Ok(Some(value));
            }
            if let IntegralKey::Legacy { kind, nums } = &identity {
                if let Some(value) = typed_alias(connection, kind, nums, typed_already_tried)? {
                    return Ok(Some(value));
                }
            }
        }
        self.old_value(connection, key)
    }

    fn descendant(
        &self,
        connection: &Connection,
        g: u32,
        d: &[u32],
        ch: &[(u32, u32)],
        typed_already_tried: Option<&IntegralKey>,
    ) -> Result<Option<Q>, String> {
        if self.has_values {
            let n = length(d.len())?;
            let mut body = Vec::new();
            for number in [g, n, length(ch.len())?]
                .into_iter()
                .chain(d.iter().copied())
                .chain(ch.iter().flat_map(|&(m, p)| [m, p]))
            {
                put(u64::from(number), &mut body);
            }
            let mut key = DESCENDANT.to_vec();
            // These fields were u32 arithmetic in the historical release
            // writer. The full identity follows them, so wrapped ranks do not
            // merge different mathematical states.
            let weight = ch
                .iter()
                .fold(0u32, |sum, &(m, p)| sum.wrapping_add(m.wrapping_mul(p)));
            key.extend_from_slice(&weight.to_be_bytes());
            key.extend_from_slice(
                &g.wrapping_mul(2)
                    .wrapping_add(n)
                    .wrapping_sub(2)
                    .to_be_bytes(),
            );
            key.extend_from_slice(&body);
            if let Some(value) = self.blob(connection, &key)? {
                return Ok(Some(value));
            }
            key.clear();
            key.extend_from_slice(OLD_DESCENDANT);
            key.extend_from_slice(&body);
            if let Some(value) = self.blob(connection, &key)? {
                return Ok(Some(value));
            }
        }
        if ch.is_empty() && stable_top(g, d.len(), d.iter().map(|&power| u64::from(power)).sum()) {
            let mut powers = d.to_vec();
            powers.sort_unstable();
            return self.legacy_get_unless_tried(
                connection,
                &old_key("st", &powers),
                typed_already_tried,
            );
        }
        Ok(None)
    }

    #[cfg(test)]
    pub(super) fn get(
        &mut self,
        connection: &mut Connection,
        key: &IntegralKey,
    ) -> Result<Option<Q>, String> {
        self.get_inner(connection, key, false)
    }

    /// The ordinary ranked/current point read has already missed. Descendants
    /// still test tag 1 before their exact aliases are marked fully examined.
    pub(super) fn get_after_current_miss(
        &mut self,
        connection: &mut Connection,
        key: &IntegralKey,
    ) -> Result<Option<Q>, String> {
        self.get_inner(connection, key, true)
    }

    fn get_inner(
        &mut self,
        connection: &mut Connection,
        key: &IntegralKey,
        current_already_tried: bool,
    ) -> Result<Option<Q>, String> {
        let typed_already_tried = current_already_tried.then_some(key);
        match key {
            IntegralKey::Legacy { kind, nums } => {
                if let Some(value) =
                    current_value_unless_tried(connection, key, typed_already_tried)?
                {
                    return Ok(Some(value));
                }
                if let Some(value) = typed_alias(connection, kind, nums, typed_already_tried)? {
                    return Ok(Some(value));
                }
                self.old_value(connection, &old_key(kind, nums))
            }
            IntegralKey::Descendant { g, d, ch } => {
                if let Some(value) = unranked_value(connection, key)? {
                    return Ok(Some(value));
                }
                self.descendant(connection, *g, d, ch, typed_already_tried)
            }
            IntegralKey::Integral {
                g,
                n,
                kappa,
                psi,
                lambda,
                ch,
            } => {
                let candidate = klp_key(*g, *n, kappa, psi, lambda, ch)?;
                if let Some(value) = self.old_value(connection, &candidate)? {
                    return Ok(Some(value));
                }
                if lambda.is_empty() {
                    let mut nums = vec![
                        *g,
                        *n,
                        pair_numbers(kappa.len())?,
                        pair_numbers(psi.len())?,
                        pair_numbers(ch.len())?,
                    ];
                    nums.extend(kappa.iter().chain(psi).chain(ch).flat_map(|&(i, p)| [i, p]));
                    if let Some(value) = self.legacy_get_unless_tried(
                        connection,
                        &old_key("ic", &nums),
                        typed_already_tried,
                    )? {
                        return Ok(Some(value));
                    }
                    if kappa.is_empty() && ch.is_empty() {
                        if let Some(d) = pure_root(*g, *n, psi)? {
                            if let Some(value) =
                                self.descendant(connection, *g, &d, &[], typed_already_tried)?
                            {
                                return Ok(Some(value));
                            }
                        }
                    }
                }
                if self.has_values && !self.klp_imported {
                    self.import_klp(connection)?;
                    return current_value(connection, key);
                }
                Ok(None)
            }
        }
    }

    fn import_klp(&mut self, connection: &mut Connection) -> Result<(), String> {
        let mut cursor: Option<Vec<u8>> = None;
        let mut encoded_key = Vec::new();
        let mut encoded_value = Vec::new();
        loop {
            // End the read statement and cursor before starting a write TX.
            let rows = {
                let mut statement = connection
                    .prepare_cached(
                        "SELECT k,v FROM faber_values WHERE k>=?1 AND k<?2
                     AND (?3 IS NULL OR k>?3) ORDER BY k LIMIT ?4",
                    )
                    .map_err(sql_error)?;
                let mapped = statement
                    .query_map(
                        params![KLP, KLP_END, cursor.as_deref(), IMPORT_ROWS as i64],
                        |row| Ok((row.get::<_, Vec<u8>>(0)?, row.get::<_, Vec<u8>>(1)?)),
                    )
                    .map_err(sql_error)?;
                mapped
                    .collect::<rusqlite::Result<Vec<_>>>()
                    .map_err(sql_error)?
            };
            if rows.is_empty() {
                self.klp_imported = true;
                return Ok(());
            }
            debug_assert!(rows.len() <= IMPORT_ROWS);
            let last = rows.last().unwrap().0.clone();
            let transaction = connection.transaction().map_err(sql_error)?;
            {
                let mut insert = transaction
                    .prepare_cached(
                        "INSERT INTO faber_integrals(k,v) VALUES(?1,?2) ON CONFLICT(k) DO NOTHING",
                    )
                    .map_err(sql_error)?;
                for (raw_key, bytes) in rows {
                    let key = canonical_klp(&raw_key)?;
                    let value = old_rational(&bytes)?;
                    encode_key_into(&key, &mut encoded_key);
                    encoded_value.clear();
                    encode_rational_into(&value, &mut encoded_value);
                    insert
                        .execute(params![encoded_key, encoded_value])
                        .map_err(sql_error)?;
                }
            }
            transaction.commit().map_err(sql_error)?;
            cursor = Some(last);
        }
    }
}

fn pair_numbers(count: usize) -> Result<u32, String> {
    length(count.checked_mul(2).ok_or("legacy pair count overflow")?)
}

fn klp_key(
    g: u32,
    n: u32,
    kappa: &[(u32, u32)],
    psi: &[(u32, u32)],
    lambda: &[(u32, u32)],
    ch: &[(u32, u32)],
) -> Result<Vec<u8>, String> {
    let mut nums = vec![
        g,
        n,
        length(kappa.len())?,
        length(psi.len())?,
        length(lambda.len())?,
        length(ch.len())?,
    ];
    nums.extend(
        kappa
            .iter()
            .chain(psi)
            .chain(lambda)
            .chain(ch)
            .flat_map(|&(i, p)| [i, p]),
    );
    Ok(old_key("klp", &nums))
}

fn stable_top(g: u32, n: usize, degree: u64) -> bool {
    let n = n as i128;
    2 * i128::from(g) - 2 + n > 0 && i128::from(degree) == 3 * i128::from(g) - 3 + n
}

fn pure_root(g: u32, n: u32, psi: &[(u32, u32)]) -> Result<Option<Vec<u32>>, String> {
    let degree: u64 = psi.iter().map(|&(_, power)| u64::from(power)).sum();
    if !stable_top(g, n as usize, degree) {
        return Ok(None);
    }
    let mut d = Vec::new();
    d.try_reserve_exact(n as usize)
        .map_err(|error| format!("legacy pure psi identity: {error}"))?;
    d.resize(n as usize, 0u32);
    for &(leg, power) in psi {
        if leg == 0 || leg > n {
            return Err("invalid legacy psi marking".into());
        }
        d[leg as usize - 1] = d[leg as usize - 1]
            .checked_add(power)
            .ok_or("legacy psi exponent overflow")?;
    }
    d.sort_unstable();
    Ok(Some(d))
}

fn normalized(mut pairs: Vec<(u32, u32)>) -> Result<Vec<(u32, u32)>, String> {
    pairs.retain(|&(_, power)| power != 0);
    pairs.sort_unstable_by_key(|&(index, _)| index);
    let mut result: Vec<(u32, u32)> = Vec::with_capacity(pairs.len());
    for (index, power) in pairs {
        if let Some(last) = result.last_mut().filter(|last| last.0 == index) {
            last.1 = last
                .1
                .checked_add(power)
                .ok_or("legacy exponent sum exceeds u32")?;
        } else {
            result.push((index, power));
        }
    }
    Ok(result)
}

fn canonical_klp(raw: &[u8]) -> Result<IntegralKey, String> {
    let nums = read_numbers(raw.strip_prefix(KLP).ok_or("invalid legacy KLP prefix")?)?;
    if nums.len() < 6 {
        return Err("truncated legacy KLP identity".into());
    }
    let (g, n) = (nums[0], nums[1]);
    if 2 * i64::from(g) - 2 + i64::from(n) <= 0 {
        return Err("unstable legacy KLP identity".into());
    }
    let count = nums[2..6]
        .iter()
        .try_fold(0usize, |total, &count| total.checked_add(count as usize))
        .and_then(|count| count.checked_mul(2))
        .and_then(|count| count.checked_add(6))
        .ok_or("legacy KLP component lengths overflow")?;
    if count != nums.len() {
        return Err("invalid legacy KLP component lengths".into());
    }
    let mut offset = 6;
    let mut component = |count: u32| {
        let end = offset + 2 * count as usize;
        let pairs = nums[offset..end]
            .chunks_exact(2)
            .map(|pair| (pair[0], pair[1]))
            .collect();
        offset = end;
        pairs
    };
    let kappa = component(nums[2]);
    let psi: Vec<(u32, u32)> = component(nums[3]);
    let lambda = component(nums[4]);
    let ch = component(nums[5]);
    if psi.iter().any(|&(leg, _)| leg == 0 || leg > n) {
        return Err("invalid legacy KLP psi marking".into());
    }
    Ok(IntegralKey::Integral {
        g,
        n,
        kappa: normalized(kappa)?,
        psi: normalized(psi)?,
        lambda: normalized(lambda)?,
        ch: normalized(ch)?,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use sc_num::{One, Zero};

    fn database(old: bool) -> Connection {
        let connection = Connection::open_in_memory().unwrap();
        connection
            .execute_batch(
                "CREATE TABLE st(k BLOB PRIMARY KEY,v TEXT NOT NULL) WITHOUT ROWID;
             CREATE TABLE faber_integrals(k BLOB PRIMARY KEY,v BLOB NOT NULL) WITHOUT ROWID;",
            )
            .unwrap();
        if old {
            connection
                .execute_batch(
                    "CREATE TABLE faber_values(k BLOB PRIMARY KEY,v BLOB NOT NULL) WITHOUT ROWID",
                )
                .unwrap();
        }
        connection
    }

    fn blob(value: &Q) -> Vec<u8> {
        let mut bytes = vec![u8::from(value < &Q::ZERO), 0, 0, 0, 0];
        for limb in value.numerator_ref().limbs() {
            bytes.extend_from_slice(&limb.to_le_bytes());
        }
        while bytes.len() > 5 && bytes.last() == Some(&0) {
            bytes.pop();
        }
        let length = u32::try_from(bytes.len() - 5).unwrap();
        bytes[1..5].copy_from_slice(&length.to_le_bytes());
        for limb in value.denominator_ref().limbs() {
            bytes.extend_from_slice(&limb.to_le_bytes());
        }
        while bytes.last() == Some(&0) {
            bytes.pop();
        }
        bytes
    }

    fn put_old(connection: &Connection, key: &[u8], value: &Q) {
        connection
            .execute(
                "INSERT INTO faber_values VALUES(?1,?2)",
                params![key, blob(value)],
            )
            .unwrap();
    }

    fn put_current(connection: &Connection, key: &IntegralKey, value: &Q) {
        let mut k = Vec::new();
        encode_key_into(key, &mut k);
        let mut v = Vec::new();
        encode_rational_into(value, &mut v);
        connection
            .execute("INSERT INTO faber_integrals VALUES(?1,?2)", params![k, v])
            .unwrap();
    }

    fn put_unranked(connection: &Connection, key: &IntegralKey, value: &Q) {
        let mut k = Vec::new();
        encode_unranked_key_into(key, &mut k);
        let mut v = Vec::new();
        encode_rational_into(value, &mut v);
        connection
            .execute("INSERT INTO faber_integrals VALUES(?1,?2)", params![k, v])
            .unwrap();
    }

    fn pure(g: u32, n: u32, psi: Vec<(u32, u32)>) -> IntegralKey {
        IntegralKey::Integral {
            g,
            n,
            kappa: vec![],
            psi,
            lambda: vec![],
            ch: vec![],
        }
    }

    #[test]
    fn old_binary_values_convert_exactly_without_decimal_roundtrips() {
        let mut large = Q::ONE;
        for _ in 0..193 {
            large *= Q::from(2);
        }
        for value in [
            Q::ZERO,
            Q::ONE,
            -Q::ONE,
            Q::from_signeds(-17, 31),
            large / Q::from(173),
        ] {
            assert_eq!(old_rational(&blob(&value)).unwrap(), value);
        }
        for bytes in [
            vec![],
            vec![2, 0, 0, 0, 0, 1],
            vec![0, 9, 0, 0, 0, 1],
            vec![0, 0, 0, 0, 0],
            vec![0, 0, 0, 0, 0, 0],
        ] {
            assert!(old_rational(&bytes).is_err());
        }
    }

    #[test]
    fn old_typed_descendant_point_fallback_keeps_c_values_and_complete_markings() {
        let mut connection = database(false);
        let key = IntegralKey::Descendant {
            g: 1,
            d: vec![0],
            ch: vec![(1, 1)],
        };
        let value = Q::from_signeds(1, 2);
        let mut old = Vec::new();
        encode_unranked_key_into(&key, &mut old);
        assert_eq!(old, [1, 1, 1, 1, 0, 1, 1, 1]);
        put_unranked(&connection, &key, &value);
        let mut reader = LegacyReader::new(&connection).unwrap();
        assert_eq!(
            reader.get(&mut connection, &key).unwrap(),
            Some(value.clone())
        );
        assert_eq!(
            reader
                .get_after_current_miss(&mut connection, &key)
                .unwrap(),
            Some(value.clone())
        );
        assert_eq!(current_value(&connection, &key).unwrap(), Some(value));
        assert_eq!(
            reader
                .get(
                    &mut connection,
                    &IntegralKey::Descendant {
                        g: 1,
                        d: vec![0, 0],
                        ch: vec![(1, 1)],
                    }
                )
                .unwrap(),
            None
        );
        assert_eq!(
            reader
                .get(
                    &mut connection,
                    &IntegralKey::Integral {
                        g: 1,
                        n: 1,
                        kappa: vec![],
                        psi: vec![],
                        lambda: vec![],
                        ch: vec![(1, 1)],
                    }
                )
                .unwrap(),
            None
        );
        let count: i64 = connection
            .query_row("SELECT count(*) FROM faber_integrals", [], |row| row.get(0))
            .unwrap();
        assert_eq!(count, 1);
        // Current data has priority; the old record is retained, not rewritten.
        put_current(&connection, &key, &Q::from_signeds(3, 7));
        assert_eq!(
            current_value(&connection, &key).unwrap(),
            Some(Q::from_signeds(3, 7))
        );
        assert_eq!(
            unranked_value(&connection, &key).unwrap(),
            Some(Q::from_signeds(1, 2))
        );
    }

    #[test]
    fn post_current_miss_preserves_zero_and_distinct_typed_alias_priority() {
        let mut connection = database(true);
        let descendant = IntegralKey::Descendant {
            g: 2,
            d: vec![2, 2, 2],
            ch: vec![],
        };
        let legacy = IntegralKey::Legacy {
            kind: "st".into(),
            nums: vec![2, 2, 2],
        };
        let raw = old_key("st", &[2, 2, 2]);
        connection
            .execute("INSERT INTO st VALUES(?1,'1/8')", [&raw])
            .unwrap();
        put_old(&connection, &raw, &Q::from_signeds(1, 12));
        put_current(&connection, &legacy, &Q::from_signeds(7, 240));
        put_unranked(&connection, &descendant, &Q::ZERO);
        let mut reader = LegacyReader::new(&connection).unwrap();
        // Tag 1 precedes every older format and a valid zero is a hit.
        assert_eq!(
            reader
                .get_after_current_miss(&mut connection, &descendant)
                .unwrap(),
            Some(Q::ZERO)
        );

        let mut unranked = Vec::new();
        encode_unranked_key_into(&descendant, &mut unranked);
        connection
            .execute("DELETE FROM faber_integrals WHERE k=?1", [&unranked])
            .unwrap();
        // Skipping the already missed descendant must still read the distinct
        // modern ST identity before its old binary/text representations.
        assert_eq!(
            reader
                .get_after_current_miss(&mut connection, &descendant)
                .unwrap(),
            Some(Q::from_signeds(7, 240))
        );
        connection
            .execute("DELETE FROM faber_integrals", [])
            .unwrap();
        put_current(&connection, &descendant, &Q::from_signeds(3, 17));
        // Conversely a Legacy request has not tried its typed descendant alias.
        assert_eq!(
            reader
                .get_after_current_miss(&mut connection, &legacy)
                .unwrap(),
            Some(Q::from_signeds(3, 17))
        );
        connection
            .execute("UPDATE faber_integrals SET v=x'00'", [])
            .unwrap();
        assert!(
            reader
                .get_after_current_miss(&mut connection, &legacy)
                .is_err()
        );
    }

    #[test]
    fn post_current_miss_keeps_corrupt_tag1_and_old_text_as_errors() {
        let mut connection = database(false);
        let descendant = IntegralKey::Descendant {
            g: 2,
            d: vec![2, 2, 2],
            ch: vec![],
        };
        put_unranked(&connection, &descendant, &Q::ONE);
        connection
            .execute("UPDATE faber_integrals SET v=x'00'", [])
            .unwrap();
        let mut reader = LegacyReader::new(&connection).unwrap();
        assert!(
            reader
                .get_after_current_miss(&mut connection, &descendant)
                .is_err()
        );
        connection
            .execute("DELETE FROM faber_integrals", [])
            .unwrap();
        connection
            .execute(
                "INSERT INTO st VALUES(?1,'1/0')",
                [old_key("st", &[2, 2, 2])],
            )
            .unwrap();
        assert!(
            reader
                .get_after_current_miss(&mut connection, &descendant)
                .is_err()
        );
    }

    #[test]
    fn legacy_alias_helpers_can_read_an_old_typed_pure_descendant() {
        let mut connection = database(false);
        let key = IntegralKey::Descendant {
            g: 2,
            d: vec![2, 2, 2],
            ch: vec![],
        };
        let value = Q::from_signeds(7, 240);
        put_unranked(&connection, &key, &value);
        let mut reader = LegacyReader::new(&connection).unwrap();
        let legacy = IntegralKey::Legacy {
            kind: "st".into(),
            nums: vec![2, 2, 2],
        };
        assert_eq!(
            reader.get(&mut connection, &legacy).unwrap(),
            Some(value.clone())
        );
        assert_eq!(
            reader
                .legacy_get(&connection, &old_key("st", &[2, 2, 2]))
                .unwrap(),
            Some(value.clone())
        );
        let ic = IntegralKey::Legacy {
            kind: "ic".into(),
            nums: vec![2, 3, 0, 6, 0, 1, 2, 2, 2, 3, 2],
        };
        assert_eq!(
            reader.get(&mut connection, &ic).unwrap(),
            Some(value.clone())
        );
        assert_eq!(current_value(&connection, &key).unwrap(), Some(value));
        assert!(!reader.klp_imported);
        let count: i64 = connection
            .query_row("SELECT count(*) FROM faber_integrals", [], |row| row.get(0))
            .unwrap();
        assert_eq!(count, 1);
    }

    #[test]
    fn both_descendant_layouts_keep_c_normalization_and_zero_markings() {
        for ranked in [false, true] {
            let mut connection = database(true);
            let mut key = if ranked {
                DESCENDANT.to_vec()
            } else {
                OLD_DESCENDANT.to_vec()
            };
            if ranked {
                key.extend_from_slice(&1u32.to_be_bytes());
                key.extend_from_slice(&1u32.to_be_bytes());
            }
            key.extend_from_slice(&[1, 1, 1, 0, 1, 1]);
            put_old(&connection, &key, &Q::from_signeds(1, 2));
            let mut reader = LegacyReader::new(&connection).unwrap();
            assert_eq!(
                reader
                    .get(
                        &mut connection,
                        &IntegralKey::Descendant {
                            g: 1,
                            d: vec![0],
                            ch: vec![(1, 1)]
                        }
                    )
                    .unwrap(),
                Some(Q::from_signeds(1, 2))
            );
            let root = IntegralKey::Integral {
                g: 1,
                n: 1,
                kappa: vec![],
                psi: vec![],
                lambda: vec![],
                ch: vec![(1, 1)],
            };
            assert_eq!(reader.get(&mut connection, &root).unwrap(), None);
        }
    }

    #[test]
    fn legacy_new_binary_old_binary_and_text_have_explicit_priority() {
        let mut connection = database(true);
        let raw = old_key("st", &[2, 2, 2]);
        connection
            .execute("INSERT INTO st VALUES(?1,'7/240')", [&raw])
            .unwrap();
        put_old(&connection, &raw, &Q::from_signeds(1, 12));
        let identity = IntegralKey::Legacy {
            kind: "st".into(),
            nums: vec![2, 2, 2],
        };
        let mut reader = LegacyReader::new(&connection).unwrap();
        assert_eq!(
            reader.legacy_get(&connection, &raw).unwrap(),
            Some(Q::from_signeds(1, 12))
        );
        put_current(&connection, &identity, &Q::from_signeds(1, 6));
        assert_eq!(
            reader.get(&mut connection, &identity).unwrap(),
            Some(Q::from_signeds(1, 6))
        );
        assert_eq!(
            reader.legacy_get(&connection, &raw).unwrap(),
            Some(Q::from_signeds(1, 6))
        );
        connection
            .execute("DELETE FROM faber_integrals", [])
            .unwrap();
        connection
            .execute("UPDATE faber_values SET v=x'00'", [])
            .unwrap();
        assert!(reader.legacy_get(&connection, &raw).is_err());
        connection.execute("DELETE FROM faber_values", []).unwrap();
        assert_eq!(
            reader.legacy_get(&connection, &raw).unwrap(),
            Some(Q::from_signeds(7, 240))
        );
        connection.execute("UPDATE st SET v='1/0'", []).unwrap();
        assert!(reader.legacy_get(&connection, &raw).is_err());
    }

    #[test]
    fn stable_pure_aliases_do_not_import_the_empty_st_convention() {
        let mut connection = database(false);
        connection
            .execute(
                "INSERT INTO st VALUES(?1,'7/240')",
                [old_key("st", &[2, 2, 2])],
            )
            .unwrap();
        connection
            .execute("INSERT INTO st VALUES(?1,'1/1')", [old_key("st", &[])])
            .unwrap();
        let mut reader = LegacyReader::new(&connection).unwrap();
        assert_eq!(
            reader
                .get(
                    &mut connection,
                    &IntegralKey::Descendant {
                        g: 2,
                        d: vec![2, 2, 2],
                        ch: vec![]
                    }
                )
                .unwrap(),
            Some(Q::from_signeds(7, 240))
        );
        assert_eq!(
            reader
                .get(&mut connection, &pure(2, 3, vec![(1, 2), (2, 2), (3, 2)]))
                .unwrap(),
            Some(Q::from_signeds(7, 240))
        );
        assert_eq!(
            reader
                .get(
                    &mut connection,
                    &IntegralKey::Descendant {
                        g: 1,
                        d: vec![],
                        ch: vec![]
                    }
                )
                .unwrap(),
            None
        );
        assert_eq!(
            reader.get(&mut connection, &pure(1, 0, vec![])).unwrap(),
            None
        );
    }

    #[test]
    fn actual_ic_uses_flat_number_counts_and_pure_roots_can_read_old_c_empty() {
        let mut connection = database(true);
        let ic = old_key("ic", &[2, 1, 2, 2, 2, 1, 2, 1, 1, 1, 1]);
        put_old(&connection, &ic, &Q::from_signeds(1, 144));
        let root = IntegralKey::Integral {
            g: 2,
            n: 1,
            kappa: vec![(1, 2)],
            psi: vec![(1, 1)],
            lambda: vec![],
            ch: vec![(1, 1)],
        };
        let mut reader = LegacyReader::new(&connection).unwrap();
        assert_eq!(
            reader.get(&mut connection, &root).unwrap(),
            Some(Q::from_signeds(1, 144))
        );
        let mut old = OLD_DESCENDANT.to_vec();
        old.extend_from_slice(&[2, 2, 0, 0, 5]);
        put_old(&connection, &old, &Q::from_signeds(1, 1152));
        assert_eq!(
            reader
                .get(&mut connection, &pure(2, 2, vec![(2, 5)]))
                .unwrap(),
            Some(Q::from_signeds(1, 1152))
        );
        assert_eq!(
            reader
                .get(&mut connection, &pure(2, 3, vec![(2, 5)]))
                .unwrap(),
            None
        );
    }

    #[test]
    fn canonical_klp_hit_does_not_start_unrelated_migration() {
        let mut connection = database(true);
        put_old(
            &connection,
            &klp_key(2, 0, &[], &[], &[(1, 3)], &[]).unwrap(),
            &Q::from_signeds(1, 2880),
        );
        put_old(&connection, &[KLP, &[128]].concat(), &Q::ONE);
        let mut reader = LegacyReader::new(&connection).unwrap();
        let key = IntegralKey::Integral {
            g: 2,
            n: 0,
            kappa: vec![],
            psi: vec![],
            lambda: vec![(1, 3)],
            ch: vec![],
        };
        assert_eq!(
            reader.get(&mut connection, &key).unwrap(),
            Some(Q::from_signeds(1, 2880))
        );
        assert!(!reader.klp_imported);
    }

    #[test]
    fn raw_klp_keys_import_in_bounded_chunks_without_deleting_history() {
        let mut connection = database(true);
        for power in 0u32..513 {
            let psi = [
                (2, 0),
                (1, power.saturating_sub(1)),
                (1, u32::from(power != 0)),
            ];
            let raw = klp_key(0, 3 + power, &[(5, 0)], &psi, &[(0, 0)], &[(1, 0)]).unwrap();
            put_old(&connection, &raw, &Q::ONE);
        }
        let mut reader = LegacyReader::new(&connection).unwrap();
        let requested = pure(0, 515, vec![(1, 512)]);
        assert_eq!(
            reader.get(&mut connection, &requested).unwrap(),
            Some(Q::ONE)
        );
        assert!(reader.klp_imported);
        let old: i64 = connection
            .query_row("SELECT count(*) FROM faber_values", [], |row| row.get(0))
            .unwrap();
        let new: i64 = connection
            .query_row("SELECT count(*) FROM faber_integrals", [], |row| row.get(0))
            .unwrap();
        assert_eq!(old, 513);
        assert_eq!(new, 513);
        assert_eq!(
            current_value(&connection, &pure(0, 3, vec![])).unwrap(),
            Some(Q::ONE)
        );
    }

    #[test]
    fn malformed_imports_fail_without_marking_the_scan_complete() {
        let mut connection = database(true);
        let malformed = klp_key(0, 3, &[], &[(4, 0)], &[], &[]).unwrap();
        put_old(&connection, &malformed, &Q::ONE);
        let mut reader = LegacyReader::new(&connection).unwrap();
        assert!(reader.get(&mut connection, &pure(0, 3, vec![])).is_err());
        assert!(!reader.klp_imported);
        assert!(
            canonical_klp(&klp_key(2, 1, &[(1, u32::MAX), (1, 1)], &[], &[], &[]).unwrap())
                .is_err()
        );
        assert!(canonical_klp(&[KLP, &[128]].concat()).is_err());
    }

    #[test]
    fn failed_later_chunk_keeps_committed_chunks_and_can_retry() {
        let mut connection = database(true);
        for power in 0u32..513 {
            let raw = klp_key(0, 3 + power, &[(5, 0)], &[(1, power)], &[], &[]).unwrap();
            put_old(&connection, &raw, &Q::ONE);
        }
        // All valid keys start with genus zero. This malformed high key sorts
        // into the later chunk, after one complete chunk has committed.
        let malformed = [KLP, &[255, 255, 255, 255, 255]].concat();
        put_old(&connection, &malformed, &Q::ONE);
        let mut reader = LegacyReader::new(&connection).unwrap();
        let requested = pure(0, 515, vec![(1, 512)]);
        assert!(reader.get(&mut connection, &requested).is_err());
        assert!(!reader.klp_imported);
        let count: i64 = connection
            .query_row("SELECT count(*) FROM faber_integrals", [], |row| row.get(0))
            .unwrap();
        assert_eq!(count, IMPORT_ROWS as i64);
        // ULEB128 byte order is not numeric order. Select the last valid raw
        // key in the same order as the import, rather than assuming power 512
        // was in the rolled-back chunk. A hit in the committed chunk may
        // correctly return without restarting the scan.
        let pending_raw: Vec<u8> = connection
            .query_row(
                "SELECT k FROM faber_values WHERE k>=?1 AND k<?2 AND k<>?3
                 ORDER BY k DESC LIMIT 1",
                params![KLP, KLP_END, malformed.as_slice()],
                |row| row.get(0),
            )
            .unwrap();
        let pending = canonical_klp(&pending_raw).unwrap();
        assert_eq!(current_value(&connection, &pending).unwrap(), None);
        connection
            .execute("DELETE FROM faber_values WHERE k=?1", [malformed])
            .unwrap();
        assert_eq!(reader.get(&mut connection, &pending).unwrap(), Some(Q::ONE));
        assert!(reader.klp_imported);
        let count: i64 = connection
            .query_row("SELECT count(*) FROM faber_integrals", [], |row| row.get(0))
            .unwrap();
        assert_eq!(count, 513);
    }

    #[test]
    fn legacy_st_can_recover_a_new_typed_pure_descendant() {
        let mut connection = database(false);
        let value = Q::from_signeds(7, 240);
        let typed = IntegralKey::Descendant {
            g: 2,
            d: vec![2, 2, 2],
            ch: vec![],
        };
        put_current(&connection, &typed, &value);
        let mut reader = LegacyReader::new(&connection).unwrap();
        let legacy = IntegralKey::Legacy {
            kind: "st".into(),
            nums: vec![2, 2, 2],
        };
        assert_eq!(
            reader.get(&mut connection, &legacy).unwrap(),
            Some(value.clone())
        );
        assert_eq!(
            reader
                .legacy_get(&connection, &old_key("st", &[2, 2, 2]))
                .unwrap(),
            Some(value)
        );
        // The exact Legacy namespace retains priority over a typed alias.
        put_current(&connection, &legacy, &Q::from_signeds(1, 12));
        assert_eq!(
            reader.get(&mut connection, &legacy).unwrap(),
            Some(Q::from_signeds(1, 12))
        );
        assert_eq!(
            current_value(&connection, &typed).unwrap(),
            Some(Q::from_signeds(7, 240))
        );
    }

    #[test]
    fn legacy_st_empty_and_nonintegral_genus_conventions_do_not_alias() {
        let mut connection = database(false);
        put_current(
            &connection,
            &IntegralKey::Descendant {
                g: 1,
                d: vec![],
                ch: vec![],
            },
            &Q::from(99),
        );
        let mut reader = LegacyReader::new(&connection).unwrap();
        let empty = IntegralKey::Legacy {
            kind: "st".into(),
            nums: vec![],
        };
        assert_eq!(reader.get(&mut connection, &empty).unwrap(), None);
        assert_eq!(
            reader.legacy_get(&connection, &old_key("st", &[])).unwrap(),
            None
        );
        assert_eq!(
            reader
                .get(
                    &mut connection,
                    &IntegralKey::Legacy {
                        kind: "st".into(),
                        nums: vec![0]
                    }
                )
                .unwrap(),
            None
        );
        put_current(&connection, &empty, &Q::ONE);
        assert_eq!(reader.get(&mut connection, &empty).unwrap(), Some(Q::ONE));
    }

    #[test]
    fn legacy_ic_normalizes_its_components_before_reading_a_typed_root() {
        let mut connection = database(false);
        let root = IntegralKey::Integral {
            g: 2,
            n: 1,
            kappa: vec![(1, 2)],
            psi: vec![(1, 1)],
            lambda: vec![],
            ch: vec![(1, 1)],
        };
        let value = Q::from_signeds(1, 144);
        put_current(&connection, &root, &value);
        let nums = vec![2, 1, 4, 4, 2, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1];
        let mut reader = LegacyReader::new(&connection).unwrap();
        let legacy = IntegralKey::Legacy {
            kind: "ic".into(),
            nums: nums.clone(),
        };
        assert_eq!(
            reader.get(&mut connection, &legacy).unwrap(),
            Some(value.clone())
        );
        assert_eq!(
            reader
                .legacy_get(&connection, &old_key("ic", &nums))
                .unwrap(),
            Some(value)
        );
        for bad in [
            vec![2, 1, 1, 0, 0, 1],
            vec![2, 1, 2, 0, 0],
            vec![2, 1, 4, 0, 0, 1, u32::MAX, 1, 1],
        ] {
            assert!(
                reader
                    .get(
                        &mut connection,
                        &IntegralKey::Legacy {
                            kind: "ic".into(),
                            nums: bad
                        }
                    )
                    .is_err()
            );
        }
    }

    #[test]
    fn legacy_ic_pure_psi_alias_preserves_n_and_never_aliases_real_ch_to_c() {
        let mut connection = database(false);
        put_current(
            &connection,
            &IntegralKey::Descendant {
                g: 2,
                d: vec![0, 5],
                ch: vec![],
            },
            &Q::from_signeds(1, 1152),
        );
        put_current(
            &connection,
            &IntegralKey::Descendant {
                g: 1,
                d: vec![0],
                ch: vec![(1, 1)],
            },
            &Q::from_signeds(1, 2),
        );
        let mut reader = LegacyReader::new(&connection).unwrap();
        let pure_ic = IntegralKey::Legacy {
            kind: "ic".into(),
            nums: vec![2, 2, 0, 2, 0, 2, 5],
        };
        assert_eq!(
            reader.get(&mut connection, &pure_ic).unwrap(),
            Some(Q::from_signeds(1, 1152))
        );
        let other_n = IntegralKey::Legacy {
            kind: "ic".into(),
            nums: vec![2, 3, 0, 2, 0, 2, 5],
        };
        assert_eq!(reader.get(&mut connection, &other_n).unwrap(), None);
        let actual_ch = IntegralKey::Legacy {
            kind: "ic".into(),
            nums: vec![1, 1, 0, 0, 2, 1, 1],
        };
        assert_eq!(reader.get(&mut connection, &actual_ch).unwrap(), None);
        put_current(
            &connection,
            &IntegralKey::Integral {
                g: 1,
                n: 1,
                kappa: vec![],
                psi: vec![],
                lambda: vec![],
                ch: vec![(1, 1)],
            },
            &Q::from_signeds(1, 24),
        );
        assert_eq!(
            reader.get(&mut connection, &actual_ch).unwrap(),
            Some(Q::from_signeds(1, 24))
        );
    }
}
