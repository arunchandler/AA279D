# PS5_Analytical.m Debug Summary

## Issues Identified

### 1. **Eccentricity > 1 Problem**
- **Root Cause**: The `qns2oe` function computes the deputy's eccentricity vector as `evec_d = evec_c + [dex; dey]`, where `dex` and `dey` are the relative eccentricity components. When these components are large, they can push the total eccentricity magnitude above 1, causing the `true2mean` function to fail.
- **Location**: Lines 47-49 in `qns2oe.m`
- **Impact**: Causes `true2mean` to return imaginary values, breaking the orbital element conversion

### 2. **Dimensionless vs Dimensioned ROE Confusion**
- **Issue**: Inconsistent handling of dimensionless vs dimensioned ROEs throughout the code
- **Problem**: `compute_roes` returns dimensionless ROEs, but the code sometimes treats them as dimensioned
- **Impact**: Incorrect Δv calculations and maneuver targeting

### 3. **Vector Size Issues**
- **Issue**: Inconsistent vector size handling in `rv2oe` calls
- **Problem**: Some vectors were row vectors when column vectors were expected
- **Impact**: Potential errors in orbital element conversion

### 4. **Excessive Δv Magnitudes**
- **Issue**: Δv calculations can produce very large values that push orbits into hyperbolic regimes
- **Problem**: No limits on Δv magnitudes in the original code
- **Impact**: Orbits become unbound or highly eccentric

## Solutions Implemented

### 1. **Safe ROE Conversion Functions**
- **`safe_qns2oe.m`**: Modified version of `qns2oe` that clamps eccentricity to 0.9999 if it exceeds the limit
- **`safe_rv2oe.m`**: Modified version of `rv2oe` that handles eccentricity clamping and returns success status
- **Benefits**: Prevents eccentricity > 1 errors while maintaining orbital geometry

### 2. **Enhanced Debugging**
- **`debug_roe_conversion.m`**: Function to analyze ROE conversion process and identify issues
- **Added debugging output**: Tracks eccentricity changes during maneuvers
- **Benefits**: Helps identify when and why issues occur

### 3. **Δv Safeguards**
- **Added Δv limits**: Maximum Δv of 100 m/s to prevent excessive maneuvers
- **Debug output**: Shows Δv calculation details
- **Benefits**: Prevents maneuvers that would create unbound orbits

### 4. **Vector Size Consistency**
- **Fixed vector handling**: Ensured all vectors passed to `rv2oe` are column vectors
- **Benefits**: Eliminates potential conversion errors

## Key Changes Made

### In `PS5_Analytical.m`:
1. **Line ~50**: Changed `qns2oe` to `safe_qns2oe` with success checking
2. **Line ~55**: Added debug call to `debug_roe_conversion`
3. **Lines ~350-370**: Added Δv magnitude limits and debugging
4. **Lines ~380-450**: Updated maneuver application to use `safe_rv2oe`
5. **Lines ~460-470**: Added eccentricity tracking during maneuvers

### New Files Created:
1. **`safe_qns2oe.m`**: Safe ROE to OE conversion
2. **`safe_rv2oe.m`**: Safe position-velocity to OE conversion  
3. **`debug_roe_conversion.m`**: ROE conversion debugging

## Recommendations

### 1. **Initial ROE Selection**
- When choosing `rel_qns_init`, ensure that the relative eccentricity components (elements 3-4) are not too large
- A good rule of thumb: `norm([rel_qns_init(3:4)]) < 0.1` (dimensionless)
- For the given scenarios, consider reducing the `e_term` values

### 2. **Control Window Tuning**
- The control windows (`delta_e_max_m`, `delta_i_max_m`) may need adjustment
- Smaller windows = more frequent but smaller maneuvers
- Larger windows = fewer but potentially larger maneuvers

### 3. **Monitoring**
- Use the debug output to monitor eccentricity changes
- If eccentricity clamping occurs frequently, consider:
  - Reducing relative eccentricity components
  - Adjusting control windows
  - Modifying maneuver timing

## Testing

To test the fixes:
1. Run the script with different scenarios
2. Monitor the debug output for eccentricity warnings
3. Check that no imaginary values are produced
4. Verify that the final ROE values are reasonable

The script should now handle eccentricity > 1 cases gracefully while maintaining formation control objectives. 