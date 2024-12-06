use pyo3::prelude::*;
use itertools::Itertools;
use contest_algorithms::math::{fft::{dft_from_reals,idft_to_reals},num::CommonField};
use std::{env, process::ExitCode};
use std::fs;
use rayon::prelude::*;
use std::sync::{Arc, Mutex};
use std::sync::atomic::{AtomicI32, Ordering};
use sysinfo::System;
use divisors;
fn main() ->ExitCode {
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
    match fs::create_dir_all("./results")
    {
        Ok(_) => println!("Results folder successfully created."),
        Err(_) => println!("Failed to create results folder, run in a folder with write permissions.")
    };
    println!("{:?}",args);
    println!("Hello, world!");
    let status = sicherman(sides);
    match status
    {
        1 => ExitCode::FAILURE,
        _ => ExitCode::SUCCESS
    }
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
fn onescombos(length : usize,onesdft : Vec<Vec<CommonField>>, maxlen : usize) -> Vec<(Vec<CommonField>, Vec<CommonField>)> {
    let all_combinations: Vec<Vec<usize>> = (0..length)
        .map(|_| vec![0, 1, 2])
        .multi_cartesian_product()
        .collect();

    all_combinations
        .into_par_iter()
        .map(|combination| {
            let first = (0..length)
                .flat_map(|i| std::iter::repeat(i).take((combination[i as usize]) as usize))
                .collect::<Vec<usize>>();
            let second = (0..length)
                .flat_map(|i| std::iter::repeat(i).take((2 - combination[i as usize]) as usize))
                .collect::<Vec<usize>>();
            (fastmuli(first, &onesdft, maxlen ), fastmuli(second, &onesdft, maxlen))
        })
        .collect()
}
fn transformall(facts : &Vec<Vec<i64>>, len : usize) -> Vec<Vec<CommonField>>
{
    facts.into_iter().map(|a| {dft_from_reals(&a, len)}).collect()
}
fn fastmuli(indices : Vec<usize>, facts : &Vec<Vec<CommonField>>, len : usize) -> Vec<CommonField>
{
    let t = indices.into_iter().map(|i| {facts[i].clone()}).collect();
    fastmul(&t,len)
}
fn fastmul(facts : &Vec<Vec<CommonField>>, len : usize) -> Vec<CommonField>
{
    if facts.len() == 0
    {
        return dft_from_reals(&vec![1], len);
    }
    facts.iter()
    .skip(1)
    .fold(facts[0].clone(), |accum, row| {
        accum.into_iter()
            .zip(row.iter())
            .map(|(acc, val)| acc * val.clone())
            .collect::<Vec<CommonField>>()
    })
}
fn muldft(a : &Vec<CommonField>, b: &Vec<CommonField>) -> Vec<CommonField>
{
    a.iter().zip(b).map(|(c,d)| {c.clone()*d.clone()}).collect()
}
fn sicherman(sides: i64) -> u32 {
    let factors1 = factorise(sides as u64);
    let maxlen : usize = ((factors1.0.iter().map(|a| {a.len()}).sum::<usize>()+factors1.1.iter().map(|a| {a.len()}).sum::<usize>()-factors1.0.len()-factors1.1.len())*2+1).next_power_of_two();
    println!("{}, {}",factors1.0.len(),factors1.1.len());
    println!("{:?}",factors1.0.len()+factors1.1.len());
    let ones = factors1.0;
    let size = usize::pow(3,ones.len() as u32)*2*maxlen*size_of::<usize>();
    let mut sys = System::new_all();
    sys.refresh_all();
    let avail = sys.total_memory() - sys.used_memory();
    if 15*avail/16 < size as u64
    {
        println!("Memory capacity is not sufficient, please use low memory version");
        return 1;
    }
    let onesdft = transformall(&ones, maxlen);
    let combos = onescombos(ones.len(),onesdft,maxlen);
    println!("combos done");
    let polyfactors = factors1.1;
    let pfdft = transformall(&polyfactors, maxlen);
    let factor_length = polyfactors.len();
    let factor_sums: Vec<i64> = polyfactors.iter().map(|x| x.iter().sum()).collect();
    
    let coeffs_list = Arc::new(Mutex::new(Vec::new()));
    let result_count = Arc::new(AtomicI32::new(0));
    let all_combinations = (0..factor_length)
        .map(|_| vec![0, 1, 2])
        .multi_cartesian_product();
    all_combinations.par_bridge().for_each(|combination| {
        let first: Vec<usize> = (0..factor_length)
            .flat_map(|i| std::iter::repeat(i).take(combination[i as usize] as usize))
            .collect();
        let product: i64 = first.iter().map(|&i| factor_sums[i]).product();
        if product != sides {
            return;
        }
        let second: Vec<usize> = (0..factor_length)
            .flat_map(|i| std::iter::repeat(i).take((2 - combination[i as usize]) as usize))
            .collect();
        let ac = fastmuli(first, &pfdft,maxlen);
        let bc = fastmuli(second,&pfdft,maxlen);
        combos.par_iter()
            .for_each(|i|
        {
            let ac2 = muldft(&i.0, &ac);
            let bc2 = muldft(&i.1, &bc);
            let ac1 = idft_to_reals(&ac2, maxlen);
            let bc1 = idft_to_reals(&bc2, maxlen);
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
        });
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
    return 0;
}