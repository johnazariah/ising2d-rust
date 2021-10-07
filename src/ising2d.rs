use std::fmt;
use rand::Rng;

struct Ising2DContext<const L: usize> {
    adj : [[[(usize, usize); 4]; L]; L],
}

impl<const L: usize> fmt::Debug for Ising2DContext<L> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {

        let mut result = String::from("Ising2DContext : \n{\n");
        for y in 0..L {
            for x in 0..L {
                result.push_str(&(format!("\t({}, {}) -> [ ", x, y)));
                for p in 0..4 {
                    let (a, b) = self.adj[x][y][p];
                    result.push_str(&(format!("({}, {}) ", a, b)));
                }
                result.push_str("]\n");
            }
        }
        result.push_str("}\n");

        return write!(f, "{}", result);
    }
}

fn initialize_context<const L: usize>() -> Ising2DContext<L> {
    let mut _adj = [[[(8888, 8888); 4]; L]; L];
    for (x, y) in itertools::iproduct!(0..L, 0..L) {
        _adj[x][y][0] = (((x + 1) % L), (y));
        _adj[x][y][1] = (((x + L - 1) % L), (y));
        _adj[x][y][2] = ((x), ((y + 1) % L));
        _adj[x][y][3] = ((x), ((y + L - 1) % L));
    }

    return Ising2DContext {
        adj : _adj
    };
}

struct Ising2D<const L: usize> {
    spins : [[bool; L]; L],
    ham : i64
}

impl<const L: usize> fmt::Debug for Ising2D<L> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {

        let mut result = String::from("Ising2D : \n{\n");
        result.push_str("spins : \n\t{\n");
        for y in 0..L {
            result.push_str("\t\t[ ");
            for x in 0..L {
                result.push_str(&(format!("{} ", if self.spins[x][y] { "+" } else { "0" })));
            }
            result.push_str("]\n");
        }
        result.push_str("\t}\n");
        result.push_str(&(format!("ham : \t{}\n", self.ham)));

        return write!(f, "{}", result);
    }
}

fn interaction_energy<const L: usize, const J : i8>(context : &Ising2DContext<L>, spins : &[[bool; L]; L], (x, y) : (usize, usize), flip_spin : bool) -> i64 {
    let spin_energy = | s1 : bool, s2 : bool | {
        return match J {
            j if j > 0  => if s1 == s2 { 1i64 } else {-1i64}
            j if j < 0 => if s1 != s2 { 1i64 } else {-1i64}
            _  => 0
        }
    };
    
    let s1 = if flip_spin { !spins[x][y] } else { spins[x][y] };

    let fold_func = | h : i64, n : usize | {
        let (nx, ny) = context.adj[x][y][n];
        let s2 = spins[nx][ny];
        return h + spin_energy(s1, s2);
    };

    return (0..3).fold(0, fold_func);
}

fn compute_hamiltonian<const L: usize, const J : i8>(context : &Ising2DContext<L>, spins : &[[bool; L]; L]) -> i64 {
    let e = itertools::iproduct!(0..L, 0..L)
        .fold(0, |_h, _cell| {return _h + interaction_energy::<L, J>(context, spins, _cell, false); });

    return -e;
}

fn initialize_lattice<const L: usize>(context : &Ising2DContext<L>) -> Ising2D<L> {
    let mut spins = [[false; L]; L];

    for (x, y) in itertools::iproduct!(0..L, 0..L) {
        spins[x][y] = rand::random();
    }
    let ham = compute_hamiltonian::<L,1>(context, &spins);

    return Ising2D { spins, ham };
}

fn evolveMC<'a, const L: usize, const J : i8>(context : &Ising2DContext<L>, lattice : &mut Ising2D<L>, T : f64, cell: (usize, usize)) {
    let accept = | dE | {
        return if dE <= 0 { 
            true
        } else {
            let probability = (-(dE as f64)/T).exp();
            let random = rand::random::<f64>();
            probability >= random
        };
    };

    let before = interaction_energy::<L,1>(context, &lattice.spins, cell, false);
    let after = interaction_energy::<L,1>(context, &lattice.spins, cell, true);
    let dE = after - before;
    
    if accept(dE) {
        lattice.ham = lattice.ham + dE;
        let (x, y) = cell;
        lattice.spins[x][y] = !lattice.spins[x][y];
    }
}

pub fn solve_mc<const L: usize>(num_iterations : i64, T : f64) {
    println!("Solving Ising2D {side} x {side} ({num_iterations} iterations)", side=L, num_iterations=num_iterations);
    let context = initialize_context::<L>();
    let mut lattice = initialize_lattice(&context);

    let mut rng = rand::thread_rng();

    println!("Before: {:?}", lattice);
    for _ in 0..num_iterations {
        let x = rng.gen_range(0..L);
        let y = rng.gen_range(0..L);
        evolveMC::<L, 1>(&context, &mut lattice, T, (x, y));
    }

    //println!("{:?}", context);
    println!("After: {:?}", lattice);
}