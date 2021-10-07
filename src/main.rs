mod ising2d;

fn main() {
    let c = ising2d::solveMC::<10>(100, 1.7);    
    println!("{:?}", c);
}
