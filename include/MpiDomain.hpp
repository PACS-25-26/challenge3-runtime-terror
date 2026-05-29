#pragma once

/**
 * @brief Structure representing the local domain assigned to a specific MPI process.
 * * This struct automatically computes the 1D row-wise domain decomposition, 
 * handling load balancing when the global grid size is not perfectly 
 * divisible by the number of MPI processes. It also accounts for the 
 * ghost cells needed for boundary communication.
 */
struct MpiDomain {
    int rank;
    int size;
    int global_n; /** @brief Total size of the global grid (n x n) */

    int local_rows; /** @brief Number of actual grid rows assigned to this MPI process */
    int start_row; /** @brief Global index (0 to n-1) of the first actual row assigned */

    /** * @brief Total rows to allocate in memory. 
     * Includes the actual rows plus 2 extra rows for ghost cells (top and bottom boundaries).
     */
    int allocated_rows; 

    /** * @brief Constructor for the MpiDomain struct.
     * * Automatically computes the 1D row-wise domain decomposition, handling load balancing when the global grid size is not perfectly divisible by the number of MPI processes. It also accounts for the ghost cells needed for boundary communication.
     * @param n Total size of the global grid (n x n).
     * @param r Rank of the MPI process.
     * @param s Total number of MPI processes.
     */
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