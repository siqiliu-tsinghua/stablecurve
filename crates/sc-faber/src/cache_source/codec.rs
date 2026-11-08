//! Exact mathematical key encoding. Callers supply the canonical identity;
//! this module preserves its order, zero entries, and component boundaries.

use super::IntegralKey;

const VERSION: u64 = 1;
const DESCENDANT: u64 = 1;
const INTEGRAL: u64 = 2;
const LEGACY: u64 = 3;
const RANKED_DESCENDANT: u64 = 4;

/// Lexicographic byte order agrees with unsigned numeric order. Common ranks
/// fit one byte; larger ranks sort first by effective big-endian byte length.
fn put_order(value: u128, output: &mut Vec<u8>) {
    if value < 128 {
        output.push(value as u8);
    } else {
        let bytes = value.to_be_bytes();
        let length = (128 - value.leading_zeros()).div_ceil(8) as usize;
        output.push(0x80 + length as u8);
        output.extend_from_slice(&bytes[bytes.len() - length..]);
    }
}

fn append_descendant_prefix(g: u32, n: usize, weight: u128, output: &mut Vec<u8>) {
    // Evaluate stability before unsigned conversion: auxiliary (1,0), (0,1)
    // and (0,2) states must not become enormous positive Euler ranks.
    let euler = 2 * i128::from(g) - 2 + n as i128;
    if euler > 0 {
        put_uleb(RANKED_DESCENDANT, output);
        put_order(weight, output);
        put_order(euler as u128, output);
    } else {
        put_uleb(DESCENDANT, output);
    }
}

/// Prefix for an already canonical raw descendant identity. The caller appends
/// its original ULEB g/n/d/ch bytes unchanged; no identity vector is rebuilt.
pub(crate) fn encode_descendant_prefix_into(g: u32, n: usize, weight: u128, output: &mut Vec<u8>) {
    output.clear();
    put_uleb(VERSION, output);
    append_descendant_prefix(g, n, weight, output);
}

fn put_uleb(mut value: u64, output: &mut Vec<u8>) {
    while value >= 128 {
        output.push((value as u8 & 127) | 128);
        value >>= 7;
    }
    output.push(value as u8);
}

fn put_len(length: usize, output: &mut Vec<u8>) {
    put_uleb(
        u64::try_from(length).expect("integral key vector length exceeds u64"),
        output,
    );
}

fn put_powers(powers: &[(u32, u32)], output: &mut Vec<u8>) {
    put_len(powers.len(), output);
    for &(index, exponent) in powers {
        put_uleb(u64::from(index), output);
        put_uleb(u64::from(exponent), output);
    }
}

/// Replace the buffer contents, retaining allocation. Stable descendants use
/// tag 4 and ordered C-weight/Euler ranks before the original full identity.
/// All identity integers retain their version-1 ULEB128 layout and every vector
/// has an explicit element count. Rank is an ordering hint, not an identity.
/// Descendant `ch` denotes normalized C insertions; Integral `ch` denotes the
/// actual Chern character. Their different tags prevent semantic aliasing.
pub(crate) fn encode_key_into(key: &IntegralKey, output: &mut Vec<u8>) {
    encode(key, output, true);
}

/// Exact pre-rank typed layout, used only for point reads of historical keys.
pub(super) fn encode_unranked_key_into(key: &IntegralKey, output: &mut Vec<u8>) {
    encode(key, output, false);
}

fn encode(key: &IntegralKey, output: &mut Vec<u8>, ranked: bool) {
    output.clear();
    put_uleb(VERSION, output);
    match key {
        IntegralKey::Descendant { g, d, ch } => {
            if ranked {
                let weight = ch
                    .iter()
                    .map(|&(index, power)| u128::from(index) * u128::from(power))
                    .sum();
                append_descendant_prefix(*g, d.len(), weight, output);
            } else {
                put_uleb(DESCENDANT, output);
            }
            put_uleb(u64::from(*g), output);
            put_len(d.len(), output);
            for &degree in d {
                put_uleb(u64::from(degree), output);
            }
            put_powers(ch, output);
        }
        IntegralKey::Integral {
            g,
            n,
            kappa,
            psi,
            lambda,
            ch,
        } => {
            put_uleb(INTEGRAL, output);
            put_uleb(u64::from(*g), output);
            put_uleb(u64::from(*n), output);
            for powers in [kappa, psi, lambda, ch] {
                put_powers(powers, output);
            }
        }
        IntegralKey::Legacy { kind, nums } => {
            put_uleb(LEGACY, output);
            put_len(kind.len(), output);
            output.extend_from_slice(kind.as_bytes());
            put_len(nums.len(), output);
            for &number in nums {
                put_uleb(u64::from(number), output);
            }
        }
    }
}

/// None means an unknown/noncanonical layout and must conservatively enable
/// recovery. Some(false) proves this exact current descendant is handled before
/// Planner::node creates a DAG state; value bytes are deliberately irrelevant.
pub(super) fn descendant_needs_dag(mut bytes: &[u8]) -> Option<bool> {
    fn uleb(bytes: &mut &[u8]) -> Option<u64> {
        let mut value = 0u64;
        for shift in (0..70).step_by(7) {
            let (&byte, tail) = bytes.split_first()?;
            *bytes = tail;
            if shift == 63 && byte > 1 {
                return None;
            }
            value |= u64::from(byte & 127) << shift;
            if byte < 128 {
                return (shift == 0 || byte != 0).then_some(value);
            }
        }
        None
    }
    fn order(bytes: &mut &[u8]) -> Option<u128> {
        let (&first, tail) = bytes.split_first()?;
        *bytes = tail;
        if first < 128 {
            return Some(u128::from(first));
        }
        let length = usize::from(first - 128);
        if !(1..=16).contains(&length) || bytes.len() < length || bytes[0] == 0 {
            return None;
        }
        let mut value = 0u128;
        for &byte in &bytes[..length] {
            value = (value << 8) | u128::from(byte);
        }
        *bytes = &bytes[length..];
        (value >= 128).then_some(value)
    }
    if uleb(&mut bytes)? != VERSION || uleb(&mut bytes)? != RANKED_DESCENDANT {
        return None;
    }
    let encoded_weight = order(&mut bytes)?;
    let encoded_euler = order(&mut bytes)?;
    let g = u32::try_from(uleb(&mut bytes)?).ok()?;
    let n = usize::try_from(uleb(&mut bytes)?).ok()?;
    if n > bytes.len() {
        return None;
    }
    let mut degree = 0u128;
    let mut previous = 0u32;
    for _ in 0..n {
        let power = u32::try_from(uleb(&mut bytes)?).ok()?;
        if power < previous {
            return None;
        }
        previous = power;
        degree = degree.checked_add(u128::from(power))?;
    }
    let count = usize::try_from(uleb(&mut bytes)?).ok()?;
    if count > bytes.len() / 2 {
        return None;
    }
    let mut weight = 0u128;
    let mut last = 0u32;
    for _ in 0..count {
        let index = u32::try_from(uleb(&mut bytes)?).ok()?;
        let power = u32::try_from(uleb(&mut bytes)?).ok()?;
        if index <= last || index % 2 == 0 || power == 0 {
            return None;
        }
        last = index;
        weight = weight.checked_add(u128::from(index) * u128::from(power))?;
    }
    let euler = 2 * i128::from(g) - 2 + n as i128;
    if !bytes.is_empty() || euler <= 0 || encoded_euler != euler as u128 || encoded_weight != weight
    {
        return None;
    }
    if weight > u128::from(u32::MAX) || g.checked_mul(3).is_none() {
        return None;
    }
    let top = 3 * i128::from(g) - 3 + n as i128;
    if top < 0 || degree.checked_add(weight)? != top as u128 {
        return Some(false);
    }
    if count == 0 {
        return Some(g != 0 && n != 1);
    }
    // The same immediate genus/vanishing branches precede DAG creation in
    // Planner::node. Wide arithmetic avoids its finite u32 rank assumptions.
    let bound = (3 * i128::from(g) - 3).max(1) as u128;
    Some(!(g <= 1 || u64::from(last) >= 2 * u64::from(g) || weight > bound))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn history_classification_uses_node_branches_and_never_result_values() {
        let immediate = [
            IntegralKey::Descendant {
                g: 0,
                d: vec![0, 0, 0],
                ch: vec![],
            },
            IntegralKey::Descendant {
                g: 2,
                d: vec![4],
                ch: vec![],
            },
            IntegralKey::Descendant {
                g: 1,
                d: vec![0],
                ch: vec![(1, 1)],
            },
            IntegralKey::Descendant {
                g: 2,
                d: vec![0],
                ch: vec![(1, 4)],
            },
            IntegralKey::Descendant {
                g: 2,
                d: vec![0, 0],
                ch: vec![(5, 1)],
            },
            IntegralKey::Descendant {
                g: 3,
                d: vec![0],
                ch: vec![(1, 1)],
            },
            IntegralKey::Descendant {
                g: 1,
                d: vec![0; 128],
                ch: vec![(1, 128)],
            },
        ];
        for key in immediate {
            assert_eq!(descendant_needs_dag(&encoded(&key)), Some(false), "{key:?}");
        }
        for key in [
            IntegralKey::Descendant {
                g: 2,
                d: vec![2, 3],
                ch: vec![],
            },
            IntegralKey::Descendant {
                g: 3,
                d: vec![2, 2, 2],
                ch: vec![(1, 3)],
            },
        ] {
            assert_eq!(descendant_needs_dag(&encoded(&key)), Some(true));
        }
    }

    #[test]
    fn history_classification_treats_old_unknown_and_malformed_keys_conservatively() {
        let key = IntegralKey::Descendant {
            g: 3,
            d: vec![2, 2, 2],
            ch: vec![(1, 3)],
        };
        let bytes = encoded(&key);
        for end in 0..bytes.len() {
            assert_eq!(descendant_needs_dag(&bytes[..end]), None);
        }
        let mut old = Vec::new();
        encode_unranked_key_into(&key, &mut old);
        assert_eq!(descendant_needs_dag(&old), None);
        for malformed in [
            {
                let mut v = bytes.clone();
                v.push(0);
                v
            },
            {
                let mut v = bytes.clone();
                v[0] = 2;
                v
            },
            {
                let mut v = bytes.clone();
                v[2] = 2;
                v
            },
            {
                let mut v = bytes.clone();
                v[3] = 1;
                v
            },
            encoded(&IntegralKey::Descendant {
                g: 2,
                d: vec![4, 0],
                ch: vec![(1, 1)],
            }),
            encoded(&IntegralKey::Descendant {
                g: 2,
                d: vec![2, 2],
                ch: vec![(1, 0)],
            }),
        ] {
            assert_eq!(descendant_needs_dag(&malformed), None);
        }
        assert_eq!(descendant_needs_dag(&encoded(&integral(1))), None);
    }

    fn encoded(key: &IntegralKey) -> Vec<u8> {
        let mut bytes = Vec::new();
        encode_key_into(key, &mut bytes);
        bytes
    }

    fn integral(n: u32) -> IntegralKey {
        IntegralKey::Integral {
            g: 2,
            n,
            kappa: vec![],
            psi: vec![],
            lambda: vec![],
            ch: vec![],
        }
    }

    #[test]
    fn version_tags_and_vector_boundaries_have_exact_bytes() {
        assert_eq!(
            encoded(&IntegralKey::Descendant {
                g: 2,
                d: vec![0, 4],
                ch: vec![(1, 3), (3, 1)],
            }),
            [1, 4, 6, 4, 2, 2, 0, 4, 2, 1, 3, 3, 1]
        );
        assert_eq!(
            encoded(&IntegralKey::Integral {
                g: 3,
                n: 2,
                kappa: vec![(1, 2)],
                psi: vec![(2, 4)],
                lambda: vec![(1, 3)],
                ch: vec![(3, 1)],
            }),
            [1, 2, 3, 2, 1, 1, 2, 1, 2, 4, 1, 1, 3, 1, 3, 1]
        );
        assert_eq!(encoded(&integral(0)), [1, 2, 2, 0, 0, 0, 0, 0]);
    }

    #[test]
    fn component_lengths_distinguish_the_same_flattened_pairs() {
        let mut keys = Vec::new();
        for component in 0..4 {
            let mut key = integral(1);
            let IntegralKey::Integral {
                kappa,
                psi,
                lambda,
                ch,
                ..
            } = &mut key
            else {
                unreachable!()
            };
            let components = [kappa, psi, lambda, ch];
            components[component].push((1, 1));
            keys.push(encoded(&key));
        }
        for (i, first) in keys.iter().enumerate() {
            for second in &keys[i + 1..] {
                assert_ne!(first, second);
            }
        }
        assert_ne!(
            encoded(&IntegralKey::Descendant {
                g: 2,
                d: vec![4],
                ch: vec![(1, 1)],
            }),
            encoded(&IntegralKey::Descendant {
                g: 2,
                d: vec![4, 1, 1],
                ch: vec![],
            })
        );
    }

    #[test]
    fn zero_markings_and_explicit_n_remain_part_of_the_key() {
        assert_ne!(encoded(&integral(0)), encoded(&integral(1)));
        assert_ne!(encoded(&integral(1)), encoded(&integral(2)));
        let descendant = |d| IntegralKey::Descendant {
            g: 2,
            d,
            ch: vec![],
        };
        assert_ne!(
            encoded(&descendant(vec![4])),
            encoded(&descendant(vec![0, 4]))
        );
        assert_ne!(
            encoded(&descendant(vec![0, 4])),
            encoded(&descendant(vec![0, 0, 4]))
        );
        let mut root = integral(2);
        let IntegralKey::Integral { psi, .. } = &mut root else {
            unreachable!()
        };
        psi.push((1, 3));
        let without_zero = encoded(&root);
        let IntegralKey::Integral { psi, .. } = &mut root else {
            unreachable!()
        };
        psi.push((2, 0));
        assert_ne!(encoded(&root), without_zero);
    }

    #[test]
    fn normalized_c_and_actual_ch_or_lambda_have_distinct_namespaces() {
        let descendant = IntegralKey::Descendant {
            g: 1,
            d: vec![0],
            ch: vec![(1, 1)],
        };
        let root_ch = IntegralKey::Integral {
            g: 1,
            n: 1,
            kappa: vec![],
            psi: vec![],
            lambda: vec![],
            ch: vec![(1, 1)],
        };
        let root_lambda = IntegralKey::Integral {
            g: 1,
            n: 1,
            kappa: vec![],
            psi: vec![],
            lambda: vec![(1, 1)],
            ch: vec![],
        };
        assert_ne!(encoded(&descendant), encoded(&root_ch));
        assert_ne!(encoded(&descendant), encoded(&root_lambda));
        assert_ne!(encoded(&root_ch), encoded(&root_lambda));
    }

    #[test]
    fn encoding_neither_sorts_nor_merges_identity_fields() {
        let descendant = |d, ch| IntegralKey::Descendant { g: 2, d, ch };
        assert_ne!(
            encoded(&descendant(vec![0, 4], vec![])),
            encoded(&descendant(vec![4, 0], vec![]))
        );
        assert_ne!(
            encoded(&descendant(vec![], vec![(1, 1), (3, 1)])),
            encoded(&descendant(vec![], vec![(3, 1), (1, 1)]))
        );
        assert_ne!(
            encoded(&descendant(vec![], vec![(1, 1), (1, 2)])),
            encoded(&descendant(vec![], vec![(1, 3)]))
        );
    }

    #[test]
    fn uleb_preserves_large_u32_fields_and_vector_lengths() {
        for (value, expected) in [
            (0, vec![0]),
            (127, vec![127]),
            (128, vec![128, 1]),
            (16_383, vec![255, 127]),
            (16_384, vec![128, 128, 1]),
            (u64::from(u32::MAX), vec![255, 255, 255, 255, 15]),
            (
                u64::MAX,
                vec![255, 255, 255, 255, 255, 255, 255, 255, 255, 1],
            ),
        ] {
            let mut bytes = Vec::new();
            put_uleb(value, &mut bytes);
            assert_eq!(bytes, expected);
        }
        let maximum = IntegralKey::Descendant {
            g: u32::MAX,
            d: vec![u32::MAX],
            ch: vec![(u32::MAX, u32::MAX)],
        };
        let mut old = Vec::new();
        encode_unranked_key_into(&maximum, &mut old);
        assert_eq!(
            old,
            [
                1, 1, 255, 255, 255, 255, 15, 1, 255, 255, 255, 255, 15, 1, 255, 255, 255, 255, 15,
                255, 255, 255, 255, 15,
            ]
        );
        let mut expected = vec![
            1, 4, 0x88, 255, 255, 255, 254, 0, 0, 0, 1, 0x85, 1, 255, 255, 255, 253,
        ];
        expected.extend_from_slice(&old[2..]);
        assert_eq!(encoded(&maximum), expected);
        let bytes = encoded(&IntegralKey::Descendant {
            g: 0,
            d: vec![0; 128],
            ch: vec![],
        });
        assert_eq!(&bytes[..7], &[1, 4, 0, 126, 0, 128, 1]);
        assert_eq!(bytes.len(), 7 + 128 + 1);
        assert!(bytes[7..].iter().all(|&byte| byte == 0));
        let mut root = integral(u32::MAX);
        let IntegralKey::Integral { g, .. } = &mut root else {
            unreachable!()
        };
        *g = u32::MAX;
        assert_eq!(
            &encoded(&root)[..12],
            &[1, 2, 255, 255, 255, 255, 15, 255, 255, 255, 255, 15]
        );
    }

    #[test]
    fn reused_buffer_replaces_old_contents_and_keeps_capacity() {
        let key = integral(2);
        let expected = encoded(&key);
        let mut buffer = Vec::with_capacity(1024);
        buffer.extend_from_slice(&[255; 100]);
        let capacity = buffer.capacity();
        encode_key_into(&key, &mut buffer);
        assert_eq!(buffer, expected);
        assert_eq!(buffer.capacity(), capacity);
        encode_key_into(&key, &mut buffer);
        assert_eq!(buffer, expected);
        assert_eq!(buffer.capacity(), capacity);
    }

    #[test]
    fn ordered_ranks_preserve_numeric_order_at_byte_boundaries() {
        for (value, expected) in [
            (0, vec![0]),
            (127, vec![127]),
            (128, vec![0x81, 128]),
            (255, vec![0x81, 255]),
            (256, vec![0x82, 1, 0]),
            (65_535, vec![0x82, 255, 255]),
            (65_536, vec![0x83, 1, 0, 0]),
        ] {
            let mut actual = Vec::new();
            put_order(value, &mut actual);
            assert_eq!(actual, expected);
        }
        let values = [
            0,
            1,
            126,
            127,
            128,
            129,
            254,
            255,
            256,
            257,
            65_535,
            65_536,
            u128::from(u64::MAX),
            u128::from(u64::MAX) + 1,
            u128::MAX,
        ];
        let bytes: Vec<_> = values
            .into_iter()
            .map(|value| {
                let mut bytes = Vec::new();
                put_order(value, &mut bytes);
                bytes
            })
            .collect();
        assert!(bytes.windows(2).all(|pair| pair[0] < pair[1]));
        let keys: Vec<_> = [
            (127, 255),
            (128, 1),
            (128, 127),
            (128, 128),
            (255, 1),
            (256, 1),
        ]
        .into_iter()
        .map(|(weight, euler)| {
            encoded(&IntegralKey::Descendant {
                g: 1,
                d: vec![0; euler],
                ch: vec![(1, weight)],
            })
        })
        .collect();
        assert!(keys.windows(2).all(|pair| pair[0] < pair[1]));
    }

    #[test]
    fn u128_c_weight_and_unstable_auxiliary_keys_keep_their_full_identity() {
        let d = vec![0];
        let ch = vec![(u32::MAX, u32::MAX); 2];
        let bytes = encoded(&IntegralKey::Descendant { g: 1, d, ch });
        assert_eq!(
            &bytes[..12],
            &[1, 4, 0x89, 1, 255, 255, 255, 252, 0, 0, 0, 2]
        );
        for (g, d) in [(0, vec![]), (0, vec![0]), (0, vec![0, 0]), (1, vec![])] {
            let key = IntegralKey::Descendant { g, d, ch: vec![] };
            let mut old = Vec::new();
            encode_unranked_key_into(&key, &mut old);
            assert_eq!(encoded(&key), old);
            assert_eq!(old[1], 1);
        }
    }

    #[test]
    fn old_typed_layout_and_non_descendant_tags_remain_exact() {
        let key = IntegralKey::Descendant {
            g: 2,
            d: vec![0, 4],
            ch: vec![(1, 3), (3, 1)],
        };
        let mut bytes = vec![255; 100];
        encode_unranked_key_into(&key, &mut bytes);
        assert_eq!(bytes, [1, 1, 2, 2, 0, 4, 2, 1, 3, 3, 1]);
        for key in [
            integral(0),
            IntegralKey::Legacy {
                kind: "st".into(),
                nums: vec![0, 0, 3],
            },
        ] {
            encode_unranked_key_into(&key, &mut bytes);
            assert_eq!(bytes, encoded(&key));
        }
    }
}
