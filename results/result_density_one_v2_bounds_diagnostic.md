# Density-1 v_2 bounds — diagnostic

## Data summary

- Source: `data/v_seq_N8388608.parquet`
- Filter: odd starts, coprime to 3
- Qualifying trajectories: 2,796,202
- Filtered rows: see stdout
- Trajectory length: min=1, median=51, max=248, mean=53.3196

## TEST B failure mode (trajectories with mean_v ≤ log_2(3))

Total failing: 0 (0.0000% of qualifying trajectories)

## Numerical notes

- Empirical density per trajectory has resolution 1/L where L is trajectory length. For L=20 (typical short trajectory), resolution is 0.05; the empirical density of v ≥ k for k ≥ 5 (null 2^{-4} = 0.0625) is harder to compare reliably for short trajectories. Length-binned analysis controls for this.
