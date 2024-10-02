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
fn repeat_elements(vec: Vec<Vec<i64>>, n: usize) -> Vec<Vec<i64>> {
    vec.into_iter()
        .flat_map(|inner_vec| std::iter::repeat(inner_vec).take(n))
        .collect()
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
fn sicherman(sides : i64)
{
    let polyfactors = factorise(sides as u64);
    let polyfactors = repeat_elements(polyfactors,3);
    let factor_length = polyfactors.len();
    println!("{}", factor_length);
    let factor_sums: Vec<i64> = polyfactors.iter().map(|x| x.iter().sum()).collect();
    let mut coeffs_list : Vec<(Vec<i64>,Vec<i64>,Vec<i64>)> = Vec::new();
    let mut result_count = 0;
    for iterlen1 in 0..(1+factor_length/3)
    {
        println!("ilen: {}", iterlen1);
        for a in (0..(factor_length)).combinations(iterlen1)
        {
            let product: i64 = a.iter()
                .map(|&i| factor_sums[i])
                .product();
            if product != sides { continue;}
            let b_and_c : Vec<usize> = (0..factor_length).filter(|i| !a.contains(i)).collect();
            for iterlen2 in 0..(factor_length-iterlen1-1)
            {
                for b in b_and_c.clone().into_iter().combinations(iterlen2)
                {
                    let product: i64 = b.iter().map(|&i| factor_sums[i]).product();
                    if product != sides {
                        return;
                    }
                    let c: Vec<usize> = (0..factor_length).filter(|i| !a.contains(i) && !b.contains(&i)).collect();
                    let mut ac = polyfactors[a[0]].clone();
                    if a.len() > 1 {
                        for i in a.iter().skip(1) {
                            ac = convolution(&ac, &polyfactors[*i]);
                        }
                    }
                    let mut bc = polyfactors[b[0]].clone();
                    if b.len() > 1 {
                        for i in b.iter().skip(1) {
                            bc = convolution(&bc, &polyfactors[*i]);
                        }
                    }
                    let mut cc = polyfactors[c[0]].clone();
                    if c.len() > 1 {
                        for i in c.iter().skip(1) {
                            cc = convolution(&cc, &polyfactors[*i]);
                        }
                    }
                    
                    if ac.iter().min().unwrap() >= &0 && bc.iter().min().unwrap() >= &0 && cc.iter().min().unwrap() >= &0 && ac.iter().sum::<i64>() == sides && bc.iter().sum::<i64>() == sides && cc.iter().sum::<i64>() == sides {
                        let mut coeffs = vec![ac.clone(), bc.clone(), cc.clone()];
                        coeffs.sort();
                        let coeffs = (coeffs[0].clone(),coeffs[1].clone(),coeffs[2].clone());
                        if !coeffs_list.iter().any(|i| i==&coeffs) {
                            result_count += 1;
                            println!("{}: {},{},{}",result_count,a.len(),b.len(),c.len());
                            coeffs_list.push(coeffs);
                        }
                    }
                }
            }
        }
    }
    
    let mut results : Vec<(Vec<i64>,Vec<i64>,Vec<i64>)> = Vec::new();
    for (ac,bc,cc) in coeffs_list
    {
        let a4 = coeff_to_sides(ac);
        let b4 = coeff_to_sides(bc);
        let c4 = coeff_to_sides(cc);
        let mut res = vec![a4.clone(), b4.clone(), c4.clone()];
        res.sort();
        let res = (res[0].clone(),res[1].clone(),res[2].clone());
        results.push(res);
    }
    results.sort();
    let mut contents : Vec<String> = Vec::new();
    for res in results {
        let resstr = format!("{:?}", res);
        let mut reschr = resstr.chars();
        reschr.next();
        reschr.next_back();
        contents.push(reschr.as_str().to_string());
    }
    let contents_string = contents.join("\n");
    fs::write(format!("./results/sicherman-d{}-test.txt",sides),contents_string).expect("Unable to write file");
    
}
