//! Nuclear energy generation: pp-chain, CNO cycle, and triple-alpha processes.

use tara::nucleosynthesis;

fn main() {
    println!("=== Nuclear Energy Generation Rates ===\n");

    // Solar core conditions
    println!("--- Solar Core (ρ=150 g/cm³, T=15.7 MK, X=0.34, X_CNO=0.01) ---");
    let rho_solar = 150.0;
    let t_solar = 15.7e6;
    let x_solar = 0.34;
    let x_cno = 0.01;

    let eps_pp = nucleosynthesis::pp_chain_rate(rho_solar, x_solar, t_solar);
    let eps_cno = nucleosynthesis::cno_cycle_rate(rho_solar, x_solar, x_cno, t_solar);
    let eps_3a = nucleosynthesis::triple_alpha_rate(rho_solar, 0.64, t_solar);

    println!("pp-chain:     {eps_pp:>12.4} erg/g/s");
    println!("CNO cycle:    {eps_cno:>12.4} erg/g/s");
    println!("Triple-alpha: {eps_3a:>12.4e} erg/g/s");
    println!("Dominant: {:?}", nucleosynthesis::dominant_process(t_solar));

    // Massive star core
    println!("\n--- Massive Star Core (ρ=100, T=25 MK, X=0.70) ---");
    let t_massive = 25e6;
    let eps_pp_m = nucleosynthesis::pp_chain_rate(100.0, 0.70, t_massive);
    let eps_cno_m = nucleosynthesis::cno_cycle_rate(100.0, 0.70, 0.01, t_massive);
    println!("pp-chain:  {eps_pp_m:>12.4} erg/g/s");
    println!("CNO cycle: {eps_cno_m:>12.4} erg/g/s");
    println!("CNO/pp ratio: {:.1}×", eps_cno_m / eps_pp_m);
    println!(
        "Dominant: {:?}",
        nucleosynthesis::dominant_process(t_massive)
    );

    // Helium burning
    println!("\n--- Helium Burning (ρ=10⁴, T=200 MK, Y=0.98) ---");
    let eps_he = nucleosynthesis::triple_alpha_rate(1e4, 0.98, 200e6);
    println!("Triple-alpha: {eps_he:>12.4e} erg/g/s");
    println!("Dominant: {:?}", nucleosynthesis::dominant_process(200e6));

    // pp-chain branch fractions across temperature
    println!("\n=== pp-Chain Branch Fractions ===");
    println!("{:>8} {:>8} {:>8} {:>8}", "T (MK)", "ppI", "ppII", "ppIII");
    for t_mk in [5.0, 8.0, 10.0, 12.0, 14.0, 17.0, 20.0, 25.0] {
        let (ppi, ppii, ppiii) = nucleosynthesis::pp_branch_fractions(t_mk * 1e6);
        println!("{t_mk:>8.0} {ppi:>8.3} {ppii:>8.3} {ppiii:>8.3}");
    }

    // Temperature scan: pp vs CNO crossover
    println!("\n=== pp vs CNO Crossover ===");
    println!(
        "{:>8} {:>14} {:>14} {:>10}",
        "T (MK)", "pp (erg/g/s)", "CNO (erg/g/s)", "Dominant"
    );
    for t_mk in [10.0, 13.0, 15.0, 17.0, 19.0, 22.0, 25.0, 30.0] {
        let t = t_mk * 1e6;
        let pp = nucleosynthesis::pp_chain_rate(150.0, 0.70, t);
        let cno = nucleosynthesis::cno_cycle_rate(150.0, 0.70, 0.01, t);
        let dom = nucleosynthesis::dominant_process(t);
        println!("{t_mk:>8.0} {pp:>14.4e} {cno:>14.4e} {dom:>10?}");
    }
}
