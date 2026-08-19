# Score Lock Procedure

1. Freeze model version and config hash.
2. Freeze input parser and feature schema.
3. Generate scores on blinded sponsor data.
4. Store outputs with timestamp and hash.
5. Do not inspect labels until the lock is documented.
6. After lock, compare against sponsor baseline and unblinded outcomes.

Hard rule:
No tuning, no calibration, no feature changes, no threshold updates after seeing labels.

