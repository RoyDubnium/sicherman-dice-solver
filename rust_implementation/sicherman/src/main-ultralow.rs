use pyo3::prelude::*;
use itertools::Itertools;
use contest_algorithms::math::fft::convolution;
use std::env;
use std::fs;
use rayon::prelude::*;
use std::sync::{Arc, Mutex};
use std::sync::atomic::{AtomicI32, Ordering};
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

fn factorise(input : u64) -> (Vec<Vec<i64>>,Vec<Vec<i64>>) {
    let mut ones : Vec<Vec<i64>> = Vec::new();
    let mut others : Vec<Vec<i64>> = Vec::new();
    for i in divisors::get_divisors(input).iter()
    {
        if i != &1
        {
            let v = cyclotomic(i.clone());
            if v.iter().sum::<i64>() == 1
            {
                ones.push(v);
            }
            else 
            {
                others.push(v);
            }
        };
    }
    let v = cyclotomic(input);
    if v.iter().sum::<i64>() == 1
    {
        ones.push(v);
    }
    else 
    {
        others.push(v);
    }
    return (ones,others);
}
fn cyclotomic(input: u64) -> Vec<i64> {
    let mut result : Vec<i64> = vec![];
    Python::with_gil(|py| {
        let sympy = py.import("sympy").unwrap();
        let x = sympy.call_method1("symbols", ("x",)).unwrap();
        let c = sympy.call_method1("cyclotomic_poly", (input, &x)).unwrap();
        let p = sympy.call_method1("Poly", (c, &x)).unwrap();
        result = p.call_method0("all_coeffs").unwrap().extract().unwrap();
    });
    return result;
}
fn coeff_to_sides(coeffs : Vec<i64>) -> Vec<i64>
{
    let mut sides = Vec::new();
    for (i, &c) in coeffs.iter().enumerate() {
        if c > 0 {
            sides.extend(vec![i as i64 + 1; c as usize]);
        }
    }
    return sides;
}
fn multiplyi(a : &Vec<usize>, polyfactors : &Vec<Vec<i64>>) -> Vec<i64>
{
    if a.len() == 0
    {
        return vec![1];
    }
    let mut ac = polyfactors[a[0]].clone();
    if a.len() > 1 {
        for i in a.iter().skip(1) {
            ac = convolution(&ac, &polyfactors[*i]);
        }
    }
    return ac
}
fn sicherman(sides: i64) {
    let factors1 = factorise(sides as u64);
    println!("{:?}",factors1);
    println!("{:?}",factors1.0.len()+factors1.1.len());
    let ones = factors1.0;
    let oneslen = ones.len();
    println!("combos done");
    let polyfactors = factors1.1;
    let factor_length = polyfactors.len();
    
    let coeffs_list = Arc::new(Mutex::new(Vec::new()));
    let result_count = Arc::new(AtomicI32::new(0));
    let all_combinations = (0..factor_length)
        .map(|_| vec![0, 1, 2])
        .multi_cartesian_product();
    all_combinations.par_bridge().for_each(|combination| {
        let first: Vec<usize> = (0..factor_length)
            .flat_map(|i| std::iter::repeat(i).take(combination[i as usize] as usize))
            .collect();
        let second: Vec<usize> = (0..factor_length)
            .flat_map(|i| std::iter::repeat(i).take((2 - combination[i as usize]) as usize))
            .collect();
        let ac = multiplyi(&first,&polyfactors);
        let bc = multiplyi(&second,&polyfactors);
        if ac.iter().sum::<i64>() == sides && bc.iter().sum::<i64>() == sides {
            for combination in (0..oneslen)
            .map(|_| vec![0, 1, 2])
            .multi_cartesian_product()
            {
                let first = (0..oneslen)
                    .flat_map(|i| std::iter::repeat(i).take((combination[i as usize]) as usize))
                    .collect::<Vec<usize>>();
                let second = (0..oneslen)
                    .flat_map(|i| std::iter::repeat(i).take((2 - combination[i as usize]) as usize))
                    .collect::<Vec<usize>>();
                let a2 = multiplyi(&first, &ones);
                let b2 = multiplyi(&second, &ones);
                let ac1 = convolution(&a2, &ac);
                let bc1 = convolution(&b2, &bc);
                if ac1.iter().min().unwrap() >= &0 && bc1.iter().min().unwrap() >= &0
                {
                    let coeffs = if ac1 > bc1 {
                        (ac1.clone(), bc1.clone())
                    } else {
                        (bc1.clone(), ac1.clone())
                    };
                    let mut coeffs_list = coeffs_list.lock().unwrap();
                    if !coeffs_list.iter().any(|i| i == &coeffs) {
                        result_count.store(result_count.load(Ordering::Relaxed)+1, Ordering::Relaxed);
                        coeffs_list.push(coeffs);
                    }
                }
            }
            
        }
    });

    let mut contents: Vec<String> = Vec::new();
    let coeffs_list = Arc::try_unwrap(coeffs_list).unwrap().into_inner().unwrap();
    let mut results : Vec<(Vec<i64>,Vec<i64>)> = Vec::new();
    for (ac,bc) in coeffs_list
    {
        let a4 = coeff_to_sides(ac);
        let b4 = coeff_to_sides(bc);
        let res = if a4 > b4 {
            (a4, b4)
        } else {
            (b4, a4)
        };
        results.push(res);
    }
    results.sort();
    for res in results {
        let resstr = format!("{:?}", res);
        let mut reschr = resstr.chars();
        reschr.next();
        reschr.next_back();
        contents.push(reschr.as_str().to_string());
    }
    let contents_string = contents.join("\n");
    fs::write(format!("./results/sicherman-d{:03}-test.txt", sides), contents_string).expect("Unable to write file");
}