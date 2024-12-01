use pyo3::prelude::*;
use itertools::Itertools;
use contest_algorithms::math::fft::convolution;
use std::env;
use std::fs;
use divisors;
fn main() {
    let args: Vec<String> = env::args().collect();
    let mut sides = 8;
    if args.len() > 1
    {
        sides = match args[args.len()-1].to_string().parse::<i64>()
        {
            Ok(number) => number,
            Err(_) => 8
        };
    }
    println!("{:?}",args);
    println!("Hello, world!");
    sicherman(sides);
}

fn factorise(input : u64) -> Vec<Vec<i64>> {
    let mut result : Vec<Vec<i64>> = Vec::new();
    for i in divisors::get_divisors(input).iter()
    {
        if i != &1
        {
            result.push(cyclotomic(i.clone()));
        };
    }
    result.push(cyclotomic(input));
    return result;
}
fn cyclotomic(input: u64) -> Vec<i64> {
    let mut result : Vec<i64> = vec![];
    Python::with_gil(|py| {
        let sympy = py.import_bound("sympy").unwrap();
        let x = sympy.call_method1("symbols", ("x",)).unwrap();
        let c = sympy.call_method1("cyclotomic_poly", (input, &x)).unwrap();
        let p = sympy.call_method1("Poly", (c, &x)).unwrap();
        result = p.call_method0("all_coeffs").unwrap().extract().unwrap();
    });
    return result;
}
fn double_vec(vec: &mut Vec<Vec<i64>>) {
    vec.extend_from_within(..);
}

fn sicherman(sides : i64)
{
    let mut polyfactors = factorise(sides as u64);
    double_vec(&mut polyfactors);
    let factor_length = polyfactors.len();
    println!("{}",factor_length);
    let factor_sums : Vec<i64> = polyfactors.iter().map(|x| x.iter().sum()).collect();
    let mut results : Vec<(Vec<i64>,Vec<i64>)> = Vec::new();
    let mut found1 = false;
    for iterlen in (factor_length/2)..factor_length
    {
        let mut found2 = false;
        println!("ilen: {}", iterlen);
        for a in (0..(factor_length)).combinations(iterlen)
        {
            let product: i64 = a.iter()
                .map(|&i| factor_sums[i])
                .product();
            if product != sides { continue;}
            
            let b : Vec<usize> = (0..factor_length).filter(|i| !a.contains(i)).collect();
            let mut ac = polyfactors[a[0]].clone();
            if a.len() > 1
            {
                for i in a.iter().skip(1)
                {
                    ac = convolution(&ac,&polyfactors[*i]);
                }
            }
            let mut bc = polyfactors[b[0]].clone();
            if b.len() > 1
            {
                for i in b.iter().skip(1)
                {
                    bc = convolution(&bc,&polyfactors[*i]);
                }
            }
            if ac.iter().min().unwrap() >= &0 && bc.iter().min().unwrap() >= &0 && ac.iter().sum::<i64>() == sides && bc.iter().sum::<i64>() == sides
            {
                found1 = true;
                found2 = true;
                let mut a4 = Vec::new();
                let mut b4 = Vec::new();
                for (i, &c) in ac.iter().enumerate() {
                    if c > 0 {
                        a4.extend(vec![i as i64 + 1; c as usize]);
                    }
                }

                for (i, &c) in bc.iter().enumerate() {
                    if c > 0 {
                        b4.extend(vec![i as i64 + 1; c as usize]);
                    }
                }
                let result = if a4 > b4 {
                    (a4.clone(), b4.clone())
                } else {
                    (b4.clone(), a4.clone())
                };

                
                if !results.iter().any(|i| i==&result) {
                    println!("{}",a.len());
                    results.push(result);
                }
            }
        }
        if found1 && !found2
        {
            break;
        }
    }
    let mut contents : Vec<String> = Vec::new();
    results.sort();
    for res in results
    {
        let resstr = format!("{:?}",res);
        let mut reschr= resstr.chars();
        reschr.next();
        reschr.next_back();
        contents.push(reschr.as_str().to_string());
    }
    let contents_string = contents.join("\n");
    fs::write(format!("./results/sicherman-d{}-test.txt",sides),contents_string).expect("Unable to write file");
    
}
