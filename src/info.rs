use crate::amd::{Info, Status};

pub fn info(info: &Info) {
    println!("\nAMD results:");

    let n = info.n;
    let ndiv = info.n_div;
    let nmultsubs_ldl = info.n_mult_subs_ldl;
    let nmultsubs_lu = info.n_mult_subs_lu;
    let lnz = info.lnz;
    let lnzd = if n >= 0 && lnz >= 0 { n + lnz } else { -1 };

    // AMD return status.
    print!("    status:                                             ");
    if info.status == Status::OK {
        println!("OK");
    } else if info.status == Status::OutOfMemory {
        println!("out of memory");
    } else if info.status == Status::Invalid {
        println!("invalid matrix");
    } else if info.status == Status::OkButJumbled {
        println!("OK, but jumbled");
    } else {
        println!("unknown");
    }

    // Statistics about the input matrix.
    print!(
        "    n, dimension of A:                                  {}\n",
        n
    );
    print!(
        "    nz, number of nonzeros in A:                        {}\n",
        info.nz
    );
    print!(
        "    symmetry of A:                                      {}\n",
        info.symmetry
    );
    print!(
        "    number of nonzeros on diagonal:                     {}\n",
        info.nz_diag
    );
    print!(
        "    nonzeros in pattern of A+A' (excl. diagonal):       {}\n",
        info.nz_a_plus_at
    );
    print!(
        "    # dense rows/columns of A+A':                       {}\n",
        info.n_dense
    );

    // Statistics about AMD's behavior.
    // print!("    memory used, in bytes:                              {}\n", info.memory);
    print!(
        "    # of memory compactions:                            {}\n",
        info.n_cmp_a
    );

    // Statistics about the ordering quality.
    print!(
        "
    The following approximate statistics are for a subsequent
    factorization of A(P,P) + A(P,P)'.  They are slight upper
    bounds if there are no dense rows/columns in A+A', and become
    looser if dense rows/columns exist.\n\n"
    );

    print!(
        "    nonzeros in L (excluding diagonal):                 {}\n",
        lnz
    );
    print!(
        "    nonzeros in L (including diagonal):                 {}\n",
        lnzd
    );
    print!(
        "    # divide operations for LDL' or LU:                 {}\n",
        ndiv
    );
    print!(
        "    # multiply-subtract operations for LDL':            {}\n",
        nmultsubs_ldl
    );
    print!(
        "    # multiply-subtract operations for LU:              {}\n",
        nmultsubs_lu
    );
    print!(
        "    max nz. in any column of L (incl. diagonal):        {}\n",
        info.d_max
    );

    // Total flop counts for various factorizations.

    if n >= 0 && ndiv >= 0 && nmultsubs_ldl >= 0 && nmultsubs_lu >= 0 {
        print!(
            "
    chol flop count for real A, sqrt counted as 1 flop: {}
    LDL' flop count for real A:                         {}
    LDL' flop count for complex A:                      {}
    LU flop count for real A (with no pivoting):        {}
    LU flop count for complex A (with no pivoting):     {}\n\n",
            n + ndiv + 2 * nmultsubs_ldl,
            ndiv + 2 * nmultsubs_ldl,
            9 * ndiv + 8 * nmultsubs_ldl,
            ndiv + 2 * nmultsubs_lu,
            9 * ndiv + 8 * nmultsubs_lu
        );
    }
}
