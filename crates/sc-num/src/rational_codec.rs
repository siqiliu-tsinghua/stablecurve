//! Exact, platform-independent binary rational encoding.
//!
//! Version 1 is `version: u8`, `negative: u8`, followed by the numerator and
//! denominator magnitudes. Each magnitude is an unsigned LEB128 byte length
//! followed by little-endian bytes, with no high zero bytes. Zero has an empty
//! numerator and denominator `[1]`. Neither platform limb width nor allocator
//! representation is part of the format.

use crate::{Q, Zero};
use malachite_base::num::conversion::traits::PowerOf2Digits;
use malachite_base::num::logic::traits::SignificantBits;

const VERSION: u8 = 1;

fn length_size(mut value: u64) -> usize {
    let mut size = 1;
    while value >= 128 {
        value >>= 7;
        size += 1;
    }
    size
}

fn write_length(mut value: u64, out: &mut Vec<u8>) {
    while value >= 128 {
        out.push((value as u8 & 127) | 128);
        value >>= 7;
    }
    out.push(value as u8);
}

/// Append one exact rational encoding to a reusable buffer.
///
/// The existing buffer contents are preserved. Magnitudes are streamed from
/// borrowed limbs without decimal conversion or temporary big-integer copies.
pub fn encode_rational_into(value: &Q, out: &mut Vec<u8>) {
    let magnitudes = [value.numerator_ref(), value.denominator_ref()];
    let lengths = magnitudes.map(|natural| natural.significant_bits().div_ceil(8));
    let additional = lengths
        .iter()
        .try_fold(2usize, |bytes, &length| {
            bytes
                .checked_add(length_size(length))?
                .checked_add(usize::try_from(length).ok()?)
        })
        .expect("rational encoding exceeds addressable memory");
    out.reserve(additional);
    out.push(VERSION);
    out.push(u8::from(value < &Q::ZERO));
    for (natural, length) in magnitudes.into_iter().zip(lengths) {
        write_length(length, out);
        let mut remaining =
            usize::try_from(length).expect("rational magnitude exceeds addressable memory");
        for limb in natural.limbs() {
            let bytes = limb.to_le_bytes();
            let count = remaining.min(bytes.len());
            out.extend_from_slice(&bytes[..count]);
            remaining -= count;
        }
    }
}

fn read_byte(input: &mut &[u8]) -> Result<u8, String> {
    let (&byte, rest) = input
        .split_first()
        .ok_or_else(|| "truncated rational encoding".to_owned())?;
    *input = rest;
    Ok(byte)
}

fn read_length(input: &mut &[u8]) -> Result<usize, String> {
    let mut value = 0u64;
    for index in 0..10 {
        let byte = read_byte(input)?;
        if index == 9 && byte > 1 {
            return Err("rational magnitude length overflows u64".to_owned());
        }
        value |= u64::from(byte & 127) << (index * 7);
        if byte & 128 == 0 {
            return usize::try_from(value)
                .map_err(|_| "rational magnitude length exceeds addressable memory".to_owned());
        }
    }
    Err("rational magnitude length overflows u64".to_owned())
}

fn read_magnitude<'a>(input: &mut &'a [u8]) -> Result<&'a [u8], String> {
    let length = read_length(input)?;
    if length > input.len() {
        return Err("truncated rational magnitude".to_owned());
    }
    let (magnitude, rest) = input.split_at(length);
    *input = rest;
    Ok(magnitude)
}

/// Decode exactly one rational, rejecting incomplete or trailing data.
///
/// The denominator must be nonzero. Safe numeric construction also reduces
/// fractions and normalizes zero; no decimal parsing or internal-layout access
/// is involved. The decoder accepts redundant high zero magnitude bytes.
pub fn decode_rational(bytes: &[u8]) -> Result<Q, String> {
    let mut input = bytes;
    if read_byte(&mut input)? != VERSION {
        return Err("unsupported rational encoding version".to_owned());
    }
    let non_negative = match read_byte(&mut input)? {
        0 => true,
        1 => false,
        _ => return Err("invalid rational encoding sign".to_owned()),
    };
    let numerator_bytes = read_magnitude(&mut input)?;
    let denominator_bytes = read_magnitude(&mut input)?;
    if !input.is_empty() {
        return Err("trailing data after rational encoding".to_owned());
    }
    // The Rational constructor determines the Natural types, so sc-num needs
    // no new direct dependency on malachite-nz for binary import.
    let denominator =
        PowerOf2Digits::<u8>::from_power_of_2_digits_asc(8, denominator_bytes.iter().copied())
            .ok_or_else(|| "invalid rational denominator magnitude".to_owned())?;
    if denominator == 0u32 {
        return Err("zero denominator in rational encoding".to_owned());
    }
    let numerator =
        PowerOf2Digits::<u8>::from_power_of_2_digits_asc(8, numerator_bytes.iter().copied())
            .ok_or_else(|| "invalid rational numerator magnitude".to_owned())?;
    Ok(Q::from_sign_and_naturals(
        non_negative,
        numerator,
        denominator,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{One, pow_unsigned};

    fn encoded(value: &Q) -> Vec<u8> {
        let mut bytes = Vec::new();
        encode_rational_into(value, &mut bytes);
        bytes
    }

    fn power_of_two(exponent: u64) -> Q {
        pow_unsigned(&Q::from(2), exponent)
    }

    #[test]
    fn small_values_have_stable_exact_bytes() {
        for (value, bytes) in [
            (Q::ZERO, vec![1, 0, 0, 1, 1]),
            (Q::ONE, vec![1, 0, 1, 1, 1, 1]),
            (Q::from_signeds(-3, 5), vec![1, 1, 1, 3, 1, 5]),
            (Q::from(256), vec![1, 0, 2, 0, 1, 1, 1]),
            (Q::from_signeds(1, 256), vec![1, 0, 1, 1, 2, 0, 1]),
        ] {
            assert_eq!(encoded(&value), bytes);
            assert_eq!(decode_rational(&bytes).unwrap(), value);
        }
    }

    #[test]
    fn byte_order_is_independent_of_limb_boundaries() {
        let value = power_of_two(64) + Q::from(0x01_0203u32);
        let bytes = vec![1, 0, 9, 3, 2, 1, 0, 0, 0, 0, 0, 1, 1, 1];
        assert_eq!(encoded(&value), bytes);
        assert_eq!(decode_rational(&bytes).unwrap(), value);
        for exponent in [0, 1, 7, 8, 31, 32, 63, 64, 65, 127, 128, 129, 255, 256] {
            let magnitude = power_of_two(exponent) + Q::ONE;
            for value in [magnitude.clone(), -magnitude] {
                assert_eq!(decode_rational(&encoded(&value)).unwrap(), value);
            }
        }
    }

    #[test]
    fn arbitrary_precision_integer_and_denominator_round_trip() {
        let numerator = power_of_two(521) + Q::ONE;
        let denominator = power_of_two(607);
        for value in [
            numerator.clone(),
            -numerator.clone(),
            &numerator / &denominator,
            -numerator / denominator,
            Q::ONE / power_of_two(2048),
        ] {
            assert_eq!(decode_rational(&encoded(&value)).unwrap(), value);
        }
    }

    #[test]
    fn magnitude_lengths_use_unsigned_leb128() {
        let value = power_of_two(1024);
        let bytes = encoded(&value);
        assert_eq!(&bytes[..4], &[1, 0, 0x81, 0x01]);
        assert_eq!(bytes.len(), 135);
        assert_eq!(decode_rational(&bytes).unwrap(), value);
        let reciprocal = Q::ONE / power_of_two(2048);
        let bytes = encoded(&reciprocal);
        assert_eq!(&bytes[..6], &[1, 0, 1, 1, 0x81, 0x02]);
        assert_eq!(decode_rational(&bytes).unwrap(), reciprocal);
    }

    #[test]
    fn encoding_appends_to_a_reusable_buffer() {
        let first = Q::from_signeds(-7, 11);
        let second = Q::ONE / power_of_two(257);
        let mut bytes = vec![42, 99];
        encode_rational_into(&first, &mut bytes);
        let boundary = bytes.len();
        encode_rational_into(&second, &mut bytes);
        assert_eq!(&bytes[..2], &[42, 99]);
        assert_eq!(decode_rational(&bytes[2..boundary]).unwrap(), first);
        assert_eq!(decode_rational(&bytes[boundary..]).unwrap(), second);
        bytes.clear();
        encode_rational_into(&Q::ZERO, &mut bytes);
        assert_eq!(decode_rational(&bytes).unwrap(), Q::ZERO);
    }

    #[test]
    fn safe_import_reduces_fractions_and_normalizes_zero() {
        assert_eq!(
            decode_rational(&[1, 1, 1, 6, 1, 8]).unwrap(),
            Q::from_signeds(-3, 4)
        );
        assert_eq!(decode_rational(&[1, 1, 0, 1, 7]).unwrap(), Q::ZERO);
        assert_eq!(
            decode_rational(&[1, 0, 2, 3, 0, 2, 5, 0]).unwrap(),
            Q::from_signeds(3, 5)
        );
    }

    #[test]
    fn truncated_headers_lengths_and_magnitudes_are_errors() {
        let bytes = encoded(&(Q::ONE / power_of_two(1024)));
        for length in 0..bytes.len() {
            assert!(
                decode_rational(&bytes[..length]).is_err(),
                "length={length}"
            );
        }
        assert!(decode_rational(&[1, 0, 0x80]).is_err());
        assert!(decode_rational(&[1, 0, 0, 0x80]).is_err());
    }

    #[test]
    fn invalid_versions_signs_zero_denominators_and_tails_are_errors() {
        for bytes in [
            vec![0, 0, 0, 1, 1],
            vec![2, 0, 0, 1, 1],
            vec![1, 2, 0, 1, 1],
            vec![1, 255, 0, 1, 1],
            vec![1, 0, 1, 7, 0],
            vec![1, 0, 1, 7, 1, 0],
            vec![1, 0, 0, 2, 0, 0],
            vec![1, 0, 0, 1, 1, 0],
        ] {
            assert!(decode_rational(&bytes).is_err(), "bytes={bytes:?}");
        }
    }

    #[test]
    fn overflowing_and_unavailable_magnitude_lengths_are_errors() {
        for tail in [vec![0xff; 10], vec![0x80; 10]] {
            let mut bytes = vec![1, 0];
            bytes.extend(tail);
            assert!(decode_rational(&bytes).is_err());
        }
        let mut bytes = vec![1, 0];
        bytes.extend([0xff; 9]);
        bytes.push(1);
        assert!(decode_rational(&bytes).is_err());
    }
}
