mod ising2d;

fn main() {
    let c = ising2d::solve_mc::<10>(100000, 1.7);    
    println!("{:?}", c);
}
