#pragma once

struct MpiDomain {
    int rank;
    int size;
    int global_n; // Total size of the global grid (n x n)

    int local_rows; // Number of actual grid rows assigned to this MPI process
    int start_row; // Global index (0 to n-1) of the first actual row assigned
    
    // Total rows to allocate in memory (actual rows + 2 ghost cells for boundaries)
    int allocated_rows; 

    // Constructor computes the load balancing automatically
    MpiDomain(int n, int r, int s) : rank(r), size(s), global_n(n) {
        // 1. Load balancing: distribute rows among processes
        int base_rows = n / size;
        int remainder = n % size;

        // The first 'remainder' processes get one extra row
        if (rank < remainder) {
            local_rows = base_rows + 1;
            // Compute the global starting index
            start_row = rank * local_rows; 
        } else {
            local_rows = base_rows;
            start_row = remainder * (base_rows + 1) + (rank - remainder) * base_rows;
        }

        // 2. Ghost Cells setup
        // Each process allocates 2 extra rows (one top, one bottom) 
        // to exchange boundary data with adjacent processes.
        allocated_rows = local_rows + 2;
    }
};