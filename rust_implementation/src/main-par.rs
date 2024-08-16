use pyo3::prelude::*;
use itertools::Itertools;
use contest_algorithms::math::fft::convolution;
use std::env;
use std::fs;
use rayon::prelude::*;
use std::sync::{Arc, Mutex};
use std::sync::atomic::{AtomicBool, Ordering};
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

fn factorise(input : Vec<i64>) -> Vec<Vec<i64>> {
    let mut result : Vec<Vec<i64>> = Vec::new();
    Python::with_gil(|py| {
        let sympy = py.import_bound("sympy").unwrap();
        let x = sympy.call_method1("symbols", ("x",)).unwrap();
        let poly = sympy.call_method1("Poly", (input,&x)).unwrap();
        let polyexp = poly.call_method0("as_expr").unwrap();
        let factors = polyexp.call_method0("factor").unwrap();
        let factors_list = factors.getattr("args").unwrap();
        let factors_len = factors_list.len().unwrap();
        for factor_idx in 0usize..factors_len
        {
            let factor = factors_list.get_item(factor_idx).unwrap();
            let factor_poly = sympy.call_method1("Poly",(factor,)).unwrap();
            let coeffs : Vec<i64> = factor_poly.call_method0("all_coeffs").unwrap().extract().unwrap();
            result.push(coeffs);
        }
    });
    return result;
}
fn double_vec(vec: &mut Vec<Vec<i64>>) {
    vec.extend_from_within(..);
}

fn sicherman(sides: i64) {
    let polyvec = vec![1; sides as usize];
    let mut polyfactors = factorise(polyvec);
    double_vec(&mut polyfactors);
    let factor_length = polyfactors.len();
    println!("{}", factor_length);
    let factor_sums: Vec<i64> = polyfactors.iter().map(|x| x.iter().sum()).collect();
    
    let results = Arc::new(Mutex::new(Vec::new()));
    let found1 = Arc::new(AtomicBool::new(false));  // Persist across iterations

    for iterlen in (factor_length / 2)..factor_length {
        let results = Arc::clone(&results);
        let found1 = Arc::clone(&found1);

        let found2 = Arc::new(AtomicBool::new(false));  // Reset for each iteration
        println!("ilen: {}", iterlen);

        (0..factor_length).combinations(iterlen)
            .par_bridge() // Convert the iterator to a parallel iterator
            .for_each(|a| {
                let product: i64 = a.iter().map(|&i| factor_sums[i]).product();
                if product != sides {
                    return;
                }

                let b: Vec<usize> = (0..factor_length).filter(|i| !a.contains(i)).collect();
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

                if ac.iter().min().unwrap() >= &0 && bc.iter().min().unwrap() >= &0 && ac.iter().sum::<i64>() == sides && bc.iter().sum::<i64>() == sides {
                    found1.store(true, Ordering::Relaxed);
                    found2.store(true, Ordering::Relaxed);  // Mark found2 as true
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

                    let mut results = results.lock().unwrap();
                    if !results.iter().any(|i| i == &result) {
                        println!("{}", a.len());
                        results.push(result);
                    }
                }
            });

        // Break the outer loop if a valid result was found in a previous iteration but not in the current one
        if found1.load(Ordering::Relaxed) && !found2.load(Ordering::Relaxed) {
            break;
        }
    }

    let mut contents: Vec<String> = Vec::new();
    let mut results = Arc::try_unwrap(results).unwrap().into_inner().unwrap();
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