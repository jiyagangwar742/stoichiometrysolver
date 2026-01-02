# Stoichiometry Solver

A powerful Python tool for solving mass-to-mass stoichiometry problems with automatic equation balancing and step-by-step solutions.

## Features

- ✅ **Automatic Equation Balancing** - Balances chemical equations using matrix algebra
- ✅ **Mass-to-Mass Conversions** - Complete stoichiometric calculations
- ✅ **Complex Formula Support** - Handles parentheses and subscripts (e.g., Ca(OH)₂)
- ✅ **Step-by-Step Solutions** - Shows all work including molar masses, mole ratios, and conversions
- ✅ **Complete Periodic Table** - Includes atomic masses for all elements
- ✅ **Modular Architecture** - Clean, maintainable code structure
- ✅ **Interactive & Scriptable** - Use interactively or import as a module

## Installation

### Requirements
- Python 3.7 or higher
- NumPy

### Setup

1. **Install NumPy:**
```bash
pip install numpy
```

2. **Download the script:**
Save `stoichiometry_solver.py` to your working directory.

3. **You're ready to go!**

## Quick Start

### Run Example Problems

Simply execute the script to see worked examples:

```bash
python stoichiometry_solver.py
```

This will show:
- Example 1: Ca(OH)₂ + HCl → CaCl₂ + H₂O
- Example 2: C₃H₈ + O₂ → CO₂ + H₂O (combustion)
- Interactive mode for your own problems

### Interactive Mode

After running the script, you'll be prompted to enter your own problem:

```
Equation (e.g., H2 + O2 -> H2O): Fe + O2 -> Fe2O3
Given compound: Fe
Given mass (g): 100
Target compound: Fe2O3
```

The solver will:
1. Balance the equation automatically
2. Calculate molar masses
3. Perform stoichiometric conversions
4. Display the complete solution with all steps

## Usage Examples

### Example 1: Basic Usage

**Problem:** Given 10 grams of Ca(OH)₂, how many grams of CaCl₂ are produced?

**Input:**
```
Equation: Ca(OH)2 + HCl -> CaCl2 + H2O
Given compound: Ca(OH)2
Given mass: 10
Target compound: CaCl2
```

**Output:**
```
1. BALANCED EQUATION:
   Ca(OH)2 + 2HCl → CaCl2 + 2H2O

2. GIVEN:
   10.000 g of Ca(OH)2
   Find: mass of CaCl2

3. MOLAR MASSES:
   Ca(OH)2: 74.092 g/mol
   CaCl2: 110.978 g/mol

4. CALCULATION:
   Step 1: Convert Ca(OH)2 to moles
           10.000 g ÷ 74.092 g/mol = 0.134975 mol

   Step 2: Use mole ratio (1/1)
           0.134975 mol × (1/1) = 0.134975 mol

   Step 3: Convert to grams of CaCl2
           0.134975 mol × 110.978 g/mol = 14.980 g

5. FINAL ANSWER:
   14.980 g of CaCl2
```

### Example 2: Combustion Reaction

**Problem:** Given 25 grams of propane (C₃H₈), how many grams of CO₂ are produced?

```python
from stoichiometry_solver import StoichiometryCalculator

calc = StoichiometryCalculator()
result = calc.solve_problem(
    equation="C3H8 + O2 -> CO2 + H2O",
    given_compound="C3H8",
    given_mass=25.0,
    target_compound="CO2",
    show_steps=True
)

print(f"\nFinal Answer: {result:.2f} g")
```

### Example 3: Using as a Module

```python
from stoichiometry_solver import (
    StoichiometryCalculator,
    MolarMassCalculator,
    EquationParser,
    EquationBalancer
)

# Calculate molar mass only
molar_mass = MolarMassCalculator.calculate("H2SO4")
print(f"Molar mass of H2SO4: {molar_mass:.3f} g/mol")

# Balance an equation only
parsed = EquationParser.parse("Fe + O2 -> Fe2O3")
balanced = EquationBalancer.balance(parsed['reactants'], parsed['products'])
print(f"Balanced coefficients: {balanced}")

# Complete problem
calc = StoichiometryCalculator()
answer = calc.solve_problem(
    equation="N2 + H2 -> NH3",
    given_compound="N2",
    given_mass=50.0,
    target_compound="NH3",
    show_steps=False  # Set to True to see all steps
)
print(f"Answer: {answer:.3f} g")
```

## Module Architecture

The solver is organized into seven modular components:

### 1. **Periodic Table Data**
- Contains atomic masses for all elements
- Easy to update or extend

### 2. **Formula Parser** (`FormulaParser`)
- Parses chemical formulas into element dictionaries
- Handles parentheses with multipliers: Ca(OH)₂ → {Ca: 1, O: 2, H: 2}
- Supports nested structures

### 3. **Molar Mass Calculator** (`MolarMassCalculator`)
- Calculates molar masses from chemical formulas
- Uses periodic table data automatically

### 4. **Equation Parser** (`EquationParser`)
- Splits equations into reactants and products
- Handles both `->` and `→` arrow formats
- Cleans and normalizes input

### 5. **Equation Balancer** (`EquationBalancer`)
- Automatically balances chemical equations
- Uses matrix algebra (SVD method)
- Returns smallest integer coefficients

### 6. **Stoichiometry Solver** (`StoichiometrySolver`)
- Performs mass-to-mass conversions
- Applies mole ratios from balanced equations
- Returns detailed calculation steps

### 7. **Main Interface** (`StoichiometryCalculator`)
- User-friendly wrapper
- Coordinates all modules
- Formats and displays solutions

## Supported Features

### Chemical Formulas
- Simple formulas: H₂O, CO₂, NaCl
- Formulas with parentheses: Ca(OH)₂, Al₂(SO₄)₃
- Complex formulas: Fe₄[Fe(CN)₆]₃
- All elements in the periodic table

### Equation Types
- Synthesis reactions
- Decomposition reactions
- Single replacement
- Double replacement
- Combustion reactions
- Any balanced chemical equation

### Input Formats
- Unbalanced equations (automatic balancing)
- Flexible arrow notation: `->` or `→`
- Spaces around operators (optional)

## API Reference

### `StoichiometryCalculator.solve_problem()`

```python
solve_problem(equation, given_compound, given_mass, target_compound, show_steps=True)
```

**Parameters:**
- `equation` (str): Chemical equation (balanced or unbalanced)
- `given_compound` (str): Formula of the given substance
- `given_mass` (float): Mass of given substance in grams
- `target_compound` (str): Formula of target substance
- `show_steps` (bool): Whether to print step-by-step solution (default: True)

**Returns:**
- `float`: Mass of target substance in grams

### `MolarMassCalculator.calculate()`

```python
calculate(formula)
```

**Parameters:**
- `formula` (str): Chemical formula

**Returns:**
- `float`: Molar mass in g/mol

### `EquationBalancer.balance()`

```python
balance(reactants, products)
```

**Parameters:**
- `reactants` (List[str]): List of reactant formulas
- `products` (List[str]): List of product formulas

**Returns:**
- `Dict[str, int]`: Dictionary mapping compounds to coefficients

## Error Handling

The solver includes comprehensive error handling for common issues:

- **Invalid chemical formulas**: Raises `ValueError` with descriptive message
- **Unknown elements**: Identifies which element is not in periodic table
- **Unbalanceable equations**: Reports balancing failures
- **Missing compounds**: Validates that given/target compounds exist in equation

## Limitations

- **Equation balancing**: Very complex equations (20+ atoms, highly nonlinear) may fail
- **Redox reactions**: Does not handle half-reactions or oxidation states
- **Aqueous solutions**: Does not account for ions or solution chemistry
- **Limiting reagents**: Currently designed for single given → target conversions
- **Percent yield**: Does not calculate theoretical vs actual yield

## Future Enhancements

Potential features for future versions:
- [ ] Limiting reagent calculations
- [ ] Percent yield calculations
- [ ] Solution chemistry (molarity, dilution)
- [ ] Gas law integration (STP conditions)
- [ ] Multiple product calculations
- [ ] GUI interface
- [ ] Export solutions to PDF

## Contributing

This is a standalone educational tool. If you'd like to extend it:

1. The modular architecture makes it easy to add new features
2. Each class is independent and can be tested separately
3. Add new functionality by creating additional modules
4. Follow the existing documentation style

## Troubleshooting

### "Could not balance equation"
- Check that your equation is chemically valid
- Ensure all formulas are typed correctly
- Try simplifying complex equations

### "Unknown element: X"
- Verify element symbols are correctly capitalized (e.g., Ca not ca)
- Check for typos in chemical formulas

### Import errors
- Ensure NumPy is installed: `pip install numpy`
- Verify Python version is 3.7+

## License

This code is provided for educational purposes. Feel free to use, modify, and distribute.

## Examples Directory

Try these additional problems:

```python
calc = StoichiometryCalculator()

# Ammonia synthesis
calc.solve_problem("N2 + H2 -> NH3", "N2", 100, "NH3")

# Rusting of iron
calc.solve_problem("Fe + O2 -> Fe2O3", "Fe", 50, "Fe2O3")

# Photosynthesis
calc.solve_problem("CO2 + H2O -> C6H12O6 + O2", "CO2", 88, "C6H12O6")

# Neutralization
calc.solve_problem("HCl + NaOH -> NaCl + H2O", "HCl", 36.5, "NaCl")
```
**Happy calculating!** 🧪⚗️🔬