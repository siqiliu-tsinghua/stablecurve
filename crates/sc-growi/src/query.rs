//! 查询表达式解析——原版 `expr.cc` 关键字表的直译。
//!
//! 语法：赋值（N=/E=/G=/D=）、关键字（rational/elliptic、points/lines/…、
//! quadric/cubic/…、point/curve/surface/threefold（3fold）/…、P^n）、停用词
//! （in/on/to/thru/through/with/and，等价逗号）、类表达式（H、psi 的多项式，
//! 支持 + * ^ ( ) 整数；`类:n` 重复；`[类,重数]` 相对条件）。
//! 关键字/赋值项之后可直接跟下一项；类之后需要分隔符。

use crate::gw::{GwClass, GwSettings};
use crate::poly::{Poly, Var};

pub struct Query {
    pub settings: GwSettings,
    pub classes: Vec<GwClass>,
}

// —— 词法 ——

#[derive(Clone, Debug, PartialEq)]
enum Tok {
    Word(String), // 标识符或关键字
    Num(i64),     // 整数（含负号仅在 poly 内处理，这里只收无符号）
    Op(char),     // = + - * ^ ( ) [ ] , : /
    Stop,         // 停用词（分隔符）
}

fn stopword(w: &str) -> bool {
    matches!(w, "in" | "on" | "to" | "thru" | "through" | "with" | "and")
}

fn tokenize(s: &str) -> Result<Vec<Tok>, String> {
    let cs: Vec<char> = s.chars().collect();
    let mut out = Vec::new();
    let mut i = 0;
    while i < cs.len() {
        let c = cs[i];
        if c.is_whitespace() {
            i += 1;
        } else if c.is_ascii_digit() {
            let start = i;
            while i < cs.len() && cs[i].is_ascii_digit() {
                i += 1;
            }
            if i < cs.len() && cs[i].is_alphabetic() {
                // 数字开头的词（3fold、4fold…）：与原版词法一致，整体作为词
                while i < cs.len() && cs[i].is_alphanumeric() {
                    i += 1;
                }
                out.push(Tok::Word(cs[start..i].iter().collect()));
                continue;
            }
            let digits: String = cs[start..i].iter().collect();
            let n: i64 = digits.parse().map_err(|_| format!("数字 {digits} 过大"))?;
            out.push(Tok::Num(n));
        } else if c.is_alphabetic() {
            let mut w = String::new();
            while i < cs.len() && (cs[i].is_alphanumeric()) {
                w.push(cs[i]);
                i += 1;
            }
            let lw = w.to_ascii_lowercase();
            if stopword(&lw) {
                out.push(Tok::Stop);
            } else {
                out.push(Tok::Word(w));
            }
        } else if matches!(
            c,
            '=' | '+' | '-' | '*' | '^' | '(' | ')' | '[' | ']' | ',' | ':' | '/'
        ) {
            out.push(Tok::Op(c));
            i += 1;
        } else {
            return Err(format!("无法识别的字符 '{c}'"));
        }
    }
    Ok(out)
}

/// N/E/G/D、重数、重复次数的上限：不变量内部字段为 i16，超限须报错而非截断。
const LIMIT: i64 = i16::MAX as i64;
/// 单项式指数以 u8 存储（见 `poly::Mono`）：类表达式的总次数上限。
const MAX_EXP: u32 = u8::MAX as u32;

/// 多项式各单项式总次数的最大值。
fn total_degree(p: &Poly) -> u32 {
    p.terms()
        .iter()
        .map(|(_, m)| m.pairs().iter().map(|&(_, e)| e as u32).sum::<u32>())
        .max()
        .unwrap_or(0)
}

/// 记号的用户可读形式（错误信息用）。
fn show(t: &Option<Tok>) -> String {
    match t {
        None => "表达式结尾".into(),
        Some(Tok::Word(w)) => format!("'{w}'"),
        Some(Tok::Num(n)) => n.to_string(),
        Some(Tok::Op(c)) => format!("'{c}'"),
        Some(Tok::Stop) => "分隔词".into(),
    }
}

// —— 语法分析 ——

struct Parser {
    toks: Vec<Tok>,
    pos: usize,
    pub settings: GwSettings,
    pub classes: Vec<GwClass>,
}

impl Parser {
    fn peek(&self) -> Option<&Tok> {
        self.toks.get(self.pos)
    }

    fn next(&mut self) -> Option<Tok> {
        let t = self.toks.get(self.pos).cloned();
        if t.is_some() {
            self.pos += 1;
        }
        t
    }

    fn eat_op(&mut self, c: char) -> Result<(), String> {
        match self.next() {
            Some(Tok::Op(x)) if x == c => Ok(()),
            other => Err(format!("期望 '{c}'，得到 {}", show(&other))),
        }
    }

    fn expect_int(&mut self) -> Result<i64, String> {
        match self.next() {
            Some(Tok::Num(n)) => Ok(n),
            other => Err(format!("期望整数，得到 {}", show(&other))),
        }
    }

    /// 整数，且须落在 `[0, max]` 内（不变量内部字段为 i16，不允许静默截断）。
    fn expect_small(&mut self, what: &str, max: i64) -> Result<i64, String> {
        let n = self.expect_int()?;
        if n > max {
            return Err(format!("{what}={n} 超出支持范围（≤ {max}）"));
        }
        Ok(n)
    }

    /// 关键字 → 设置项。命中时自行消费词元（含 P^n 的 ^n 后缀）。
    fn keyword(&mut self, w: &str) -> Result<bool, String> {
        let lw = w.to_ascii_lowercase();
        let matched = matches!(
            lw.as_str(),
            "rational"
                | "elliptic"
                | "points"
                | "lines"
                | "conics"
                | "quadrics"
                | "cubics"
                | "quartics"
                | "quintics"
                | "sextics"
                | "septics"
                | "octics"
                | "nonics"
                | "quadric"
                | "cubic"
                | "quartic"
                | "quintic"
                | "sextic"
                | "septic"
                | "octic"
                | "nonic"
                | "point"
                | "curve"
                | "surface"
                | "threefold"
                | "3fold"
                | "fourfold"
                | "4fold"
                | "fivefold"
                | "5fold"
                | "sixfold"
                | "6fold"
                | "sevenfold"
                | "7fold"
                | "eightfold"
                | "8fold"
                | "ninefold"
                | "9fold"
                | "p"
        );
        if !matched {
            return Ok(false);
        }
        self.next(); // 消费关键字词
        match lw.as_str() {
            "rational" => self.settings.genus = 0,
            "elliptic" => self.settings.genus = 1,
            "points" => self.settings.cdeg = 0,
            "lines" => self.settings.cdeg = 1,
            "conics" | "quadrics" => self.settings.cdeg = 2,
            "cubics" => self.settings.cdeg = 3,
            "quartics" => self.settings.cdeg = 4,
            "quintics" => self.settings.cdeg = 5,
            "sextics" => self.settings.cdeg = 6,
            "septics" => self.settings.cdeg = 7,
            "octics" => self.settings.cdeg = 8,
            "nonics" => self.settings.cdeg = 9,
            "quadric" => self.settings.deg = 2,
            "cubic" => self.settings.deg = 3,
            "quartic" => self.settings.deg = 4,
            "quintic" => self.settings.deg = 5,
            "sextic" => self.settings.deg = 6,
            "septic" => self.settings.deg = 7,
            "octic" => self.settings.deg = 8,
            "nonic" => self.settings.deg = 9,
            "point" => self.settings.dim = 0,
            "curve" => self.settings.dim = 1,
            "surface" => self.settings.dim = 2,
            "threefold" | "3fold" => self.settings.dim = 3,
            "fourfold" | "4fold" => self.settings.dim = 4,
            "fivefold" | "5fold" => self.settings.dim = 5,
            "sixfold" | "6fold" => self.settings.dim = 6,
            "sevenfold" | "7fold" => self.settings.dim = 7,
            "eightfold" | "8fold" => self.settings.dim = 8,
            "ninefold" | "9fold" => self.settings.dim = 9,
            "p" => {
                // P^n：射影空间
                self.eat_op('^')?;
                let n = self.expect_small("P^n 的 n", LIMIT)?;
                self.settings.dim = n as i16;
                self.settings.deg = 1;
            }
            _ => unreachable!(),
        }
        Ok(true)
    }

    /// 一项：赋值 / 关键字 / 类
    /// 解析一项。返回该项是否为关键字/赋值（这类项之后不需要分隔符，
    /// 直译原版 expr.cc 的 "AGAIN: After a keyword, we don't need a comma"）。
    fn item(&mut self) -> Result<bool, String> {
        // 赋值：X= 值
        if let Some(Tok::Word(w)) = self.peek().cloned() {
            if matches!(w.as_str(), "N" | "n" | "E" | "e" | "G" | "g" | "D" | "d")
                && matches!(self.toks.get(self.pos + 1), Some(Tok::Op('=')))
            {
                self.next();
                self.next();
                let v = self.expect_small(&w.to_ascii_uppercase(), LIMIT)? as i16;
                match w.to_ascii_uppercase().as_str() {
                    "N" => self.settings.dim = v,
                    "E" => self.settings.deg = v,
                    "G" => self.settings.genus = v,
                    "D" => self.settings.cdeg = v,
                    _ => unreachable!(),
                }
                return Ok(true);
            }
            if self.keyword(&w)? {
                return Ok(true); // keyword() 已自行消费词元
            }
        }
        // 类：[poly,mult](:count)? | poly(:count)?
        let (poly, mult) = if matches!(self.peek(), Some(Tok::Op('['))) {
            self.eat_op('[')?;
            let p = self.poly()?;
            self.eat_op(',')?;
            let m = self.expect_small("重数", LIMIT)? as i16;
            self.eat_op(']')?;
            (p, m)
        } else {
            (self.poly()?, 0)
        };
        let count = if matches!(self.peek(), Some(Tok::Op(':'))) {
            self.next();
            self.expect_small("重复次数", LIMIT)? as u32
        } else {
            1
        };
        self.classes.push(GwClass { poly, mult, count });
        Ok(false)
    }

    /// 多项式：term ((+|-) term)*
    fn poly(&mut self) -> Result<Poly, String> {
        let mut p = self.term()?;
        loop {
            match self.peek() {
                Some(Tok::Op('+')) => {
                    self.next();
                    p = &p + &self.term()?;
                }
                Some(Tok::Op('-')) => {
                    self.next();
                    let t = self.term()?;
                    p = &p + &Poly::neg(&t);
                }
                _ => break,
            }
        }
        Ok(p)
    }

    fn term(&mut self) -> Result<Poly, String> {
        let mut p = self.factor()?;
        while matches!(self.peek(), Some(Tok::Op('*'))) {
            self.next();
            let f = self.factor()?;
            if total_degree(&p) + total_degree(&f) > MAX_EXP {
                return Err(format!("幂次过高（总次数须 ≤ {MAX_EXP}）"));
            }
            p = &p * &f;
        }
        Ok(p)
    }

    fn factor(&mut self) -> Result<Poly, String> {
        let base = self.base()?;
        if matches!(self.peek(), Some(Tok::Op('^'))) {
            self.next();
            let n = self.expect_small("幂次", MAX_EXP as i64)? as u32;
            if total_degree(&base) * n > MAX_EXP {
                return Err(format!("幂次过高（总次数须 ≤ {MAX_EXP}）"));
            }
            Ok(base.pow(n))
        } else {
            Ok(base)
        }
    }

    fn base(&mut self) -> Result<Poly, String> {
        match self.next() {
            Some(Tok::Num(n)) => Ok(Poly::from(n)),
            Some(Tok::Word(w)) => {
                let lw = w.to_ascii_lowercase();
                match lw.as_str() {
                    "h" => Ok(Poly::from(Var::H)),
                    "psi" | "c" => Ok(Poly::from(Var::Psi)),
                    // 原类：pc(k) = 第 k 号原始上同调类（E>1 超曲面专用）。
                    // 原版 var pc(int) → var('P',k)；输入语法为本实现扩展。
                    "pc" => {
                        self.eat_op('(')?;
                        let k = self.expect_int()?;
                        self.eat_op(')')?;
                        if k <= 0 || k > 1000 {
                            return Err(format!("原类编号 pc({k}) 无效"));
                        }
                        Ok(Poly::from(Var::Prim(k as u16)))
                    }
                    other => Err(format!("未知符号 '{other}'")),
                }
            }
            Some(Tok::Op('(')) => {
                let p = self.poly()?;
                self.eat_op(')')?;
                Ok(p)
            }
            other => Err(format!("期望类表达式，得到 {}", show(&other))),
        }
    }
}

pub fn parse(expr: &str) -> Result<Query, String> {
    let toks = tokenize(expr)?;
    let mut p = Parser {
        toks,
        pos: 0,
        settings: GwSettings::default(),
        classes: Vec::new(),
    };
    if p.peek().is_none() {
        return Err("空表达式".into());
    }
    // 关键字/赋值之后可以直接跟下一项；类之后需要分隔符
    let mut need_sep = !p.item()?;
    while let Some(t) = p.peek() {
        if matches!(t, Tok::Stop | Tok::Op(',')) {
            p.next();
            need_sep = false;
            if matches!(p.peek(), Some(Tok::Stop)) {
                continue; // 连续停用词
            }
            if p.peek().is_some() {
                need_sep = !p.item()?;
            }
        } else if !need_sep {
            need_sep = !p.item()?;
        } else {
            return Err(format!(
                "意外的记号 {}（项目之间需要分隔）",
                show(&Some(t.clone()))
            ));
        }
    }
    Ok(Query {
        settings: p.settings,
        classes: p.classes,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn digit_led_keywords() {
        // 原版 expr.cc 的 3fold/4fold… 关键字
        let q = parse("lines in quintic 3fold").unwrap();
        assert_eq!((q.settings.dim, q.settings.deg, q.settings.cdeg), (3, 5, 1));
        let q = parse("points in 4fold").unwrap();
        assert_eq!(q.settings.dim, 4);
        // 数字后接运算符仍是数字
        let q = parse("D=2,H^2:5").unwrap();
        assert_eq!(q.classes[0].count, 5);
    }

    #[test]
    fn out_of_range_is_an_error() {
        // 曾经 `as i16` 截断：D=65537 被当作 D=1
        assert!(parse("D=65537 in P^2 thru H^2:2").is_err());
        assert!(parse("P^40000").is_err());
        assert!(parse("H:99999999").is_err());
        assert!(parse("D=99999999999999999999999").is_err());
        // 单项式指数以 u8 存储：超限报错而非静默丢弃
        assert!(parse("N=2,psi^300").is_err());
        assert!(parse("N=2,psi^200*psi^100").is_err());
        assert!(parse("N=2,psi^255").is_ok());
    }

    #[test]
    fn readable_errors() {
        let Err(e) = parse("lines in P^2 thru H^2 foo") else {
            panic!("应当报错");
        };
        assert!(e.contains("'foo'"), "{e}");
        assert!(!e.contains('"'), "{e}");
    }
}
