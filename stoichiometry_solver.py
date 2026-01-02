"""
Stoichiometry Solver for Mass-to-Mass Problems
A modular tool for solving stoichiometry problems with automatic equation balancing
"""

import re
import numpy as np
from fractions import Fraction
from typing import Dict, List, Tuple


# ============================================================================
# MODULE 1: PERIODIC TABLE DATA
# ============================================================================

PERIODIC_TABLE = {
    'H': 1.008, 'He': 4.003, 'Li': 6.941, 'Be': 9.012, 'B': 10.81,
    'C': 12.011, 'N': 14.007, 'O': 15.999, 'F': 18.998, 'Ne': 20.180,
    'Na': 22.990, 'Mg': 24.305, 'Al': 26.982, 'Si': 28.085, 'P': 30.974,
    'S': 32.06, 'Cl': 35.45, 'Ar': 39.948, 'K': 39.098, 'Ca': 40.078,
    'Sc': 44.956, 'Ti': 47.867, 'V': 50.942, 'Cr': 51.996, 'Mn': 54.938,
    'Fe': 55.845, 'Co': 58.933, 'Ni': 58.693, 'Cu': 63.546, 'Zn': 65.38,
    'Ga': 69.723, 'Ge': 72.630, 'As': 74.922, 'Se': 78.971, 'Br': 79.904,
    'Kr': 83.798, 'Rb': 85.468, 'Sr': 87.62, 'Y': 88.906, 'Zr': 91.224,
    'Nb': 92.906, 'Mo': 95.95, 'Tc': 98.0, 'Ru': 101.07, 'Rh': 102.91,
    'Pd': 106.42, 'Ag': 107.87, 'Cd': 112.41, 'In': 114.82, 'Sn': 118.71,
    'Sb': 121.76, 'Te': 127.60, 'I': 126.90, 'Xe': 131.29, 'Cs': 132.91,
    'Ba': 137.33, 'La': 138.91, 'Ce': 140.12, 'Pr': 140.91, 'Nd': 144.24,
    'Pm': 145.0, 'Sm': 150.36, 'Eu': 151.96, 'Gd': 157.25, 'Tb': 158.93,
    'Dy': 163.29, 'Ho': 164.93, 'Er': 167.26, 'Tm': 168.93, 'Yb': 173.05,
    'Lu': 174.97, 'Hf': 178.49, 'Ta': 180.95, 'W': 183.84, 'Re': 186.21,
    'Os': 190.23, 'Ir': 192.22, 'Pt': 195.08, 'Au': 196.97, 'Hg': 200.59,
    'Tl': 204.38, 'Pb': 207.2, 'Bi': 208.98, 'Po': 209.0, 'At': 210.0,
    'Rn': 222.0, 'Fr': 223.0, 'Ra': 226.0, 'Ac': 227.0, 'Th': 232.04,
    'Pa': 231.04, 'U': 238.03, 'Np': 237.0, 'Pu': 244.0, 'Am': 243.0,
}


# ============================================================================
# MODULE 2: CHEMICAL FORMULA PARSER
# ============================================================================

class FormulaParser:
    """Parses chemical formulas into element composition dictionaries"""
    
    @staticmethod
    def parse(formula: str) -> Dict[str, int]:
        """
        Parse a chemical formula into elements and their counts
        
        Args:
            formula: Chemical formula string (e.g., "Ca(OH)2")
            
        Returns:
            Dictionary mapping element symbols to counts (e.g., {'Ca': 1, 'O': 2, 'H': 2})
        """
        formula = formula.strip()
        elements = {}
        
        # Expand parentheses first
        formula = FormulaParser._expand_parentheses(formula)
        
        # Pattern matches: Element symbol (uppercase + optional lowercase) followed by optional number
        pattern = r'([A-Z][a-z]?)(\d*)'
        matches = re.findall(pattern, formula)
        
        for element, count in matches:
            if element:  # Skip empty matches
                count = int(count) if count else 1
                elements[element] = elements.get(element, 0) + count
        
        return elements
    
    @staticmethod
    def _expand_parentheses(formula: str) -> str:
        """Expand parentheses in chemical formulas"""
        while '(' in formula:
            # Find innermost parentheses
            match = re.search(r'\(([^()]+)\)(\d*)', formula)
            if not match:
                break
            
            content = match.group(1)
            multiplier = int(match.group(2)) if match.group(2) else 1
            
            # Expand the contents
            expanded = ""
            pattern = r'([A-Z][a-z]?)(\d*)'
            for element, count in re.findall(pattern, content):
                count = int(count) if count else 1
                expanded += f"{element}{count * multiplier}"
            
            # Replace in formula
            formula = formula[:match.start()] + expanded + formula[match.end():]
        
        return formula


# ============================================================================
# MODULE 3: MOLAR MASS CALCULATOR
# ============================================================================

class MolarMassCalculator:
    """Calculates molar masses of chemical compounds"""
    
    @staticmethod
    def calculate(formula: str) -> float:
        """
        Calculate the molar mass of a compound
        
        Args:
            formula: Chemical formula string
            
        Returns:
            Molar mass in g/mol
        """
        elements = FormulaParser.parse(formula)
        molar_mass = 0.0
        
        for element, count in elements.items():
            if element not in PERIODIC_TABLE:
                raise ValueError(f"Unknown element: {element}")
            molar_mass += PERIODIC_TABLE[element] * count
        
        return molar_mass


# ============================================================================
# MODULE 4: EQUATION PARSER
# ============================================================================

class EquationParser:
    """Parses chemical equation strings"""
    
    @staticmethod
    def parse(equation: str) -> Dict[str, List[str]]:
        """
        Parse a chemical equation into reactants and products
        
        Args:
            equation: Chemical equation string (e.g., "H2 + O2 -> H2O")
            
        Returns:
            Dictionary with 'reactants' and 'products' lists
        """
        # Split on arrow (-> or →)
        parts = re.split(r'\s*-+>\s*|\s*→\s*', equation)
        if len(parts) != 2:
            raise ValueError("Equation must contain exactly one arrow (-> or →)")
        
        reactants_str, products_str = parts
        
        # Split each side on + and clean up
        reactants = [s.strip() for s in reactants_str.split('+')]
        products = [s.strip() for s in products_str.split('+')]
        
        # Remove any coefficient prefixes for now (we'll balance later)
        reactants = [EquationParser._remove_coefficient(r) for r in reactants]
        products = [EquationParser._remove_coefficient(p) for p in products]
        
        return {
            'reactants': reactants,
            'products': products
        }
    
    @staticmethod
    def _remove_coefficient(compound: str) -> str:
        """Remove leading numerical coefficient from compound"""
        return re.sub(r'^\d+\s*', '', compound).strip()


# ============================================================================
# MODULE 5: EQUATION BALANCER
# ============================================================================

class EquationBalancer:
    """Balances chemical equations using matrix algebra"""
    
    @staticmethod
    def balance(reactants: List[str], products: List[str]) -> Dict[str, int]:
        """
        Balance a chemical equation
        
        Args:
            reactants: List of reactant formulas
            products: List of product formulas
            
        Returns:
            Dictionary mapping each compound to its coefficient
        """
        all_compounds = reactants + products
        
        # Get all unique elements
        all_elements = set()
        for compound in all_compounds:
            elements = FormulaParser.parse(compound)
            all_elements.update(elements.keys())
        
        all_elements = sorted(list(all_elements))
        
        # Build matrix (rows = elements, cols = compounds)
        matrix = []
        for element in all_elements:
            row = []
            for i, compound in enumerate(all_compounds):
                elements = FormulaParser.parse(compound)
                count = elements.get(element, 0)
                # Reactants positive, products negative
                if i < len(reactants):
                    row.append(count)
                else:
                    row.append(-count)
            matrix.append(row)
        
        # Convert to numpy array
        matrix = np.array(matrix, dtype=float)
        
        # Find null space (solution where sum = 0 for each element)
        try:
            # Use SVD to find null space
            _, _, vh = np.linalg.svd(matrix)
            null_space = vh[-1, :]
            
            # Convert to positive integers
            # Make all positive
            if null_space[0] < 0:
                null_space = -null_space
            
            # Convert to fractions to find LCM
            fractions = [Fraction(x).limit_denominator(1000) for x in null_space]
            
            # Find LCM of denominators
            denominators = [f.denominator for f in fractions]
            lcm = denominators[0]
            for d in denominators[1:]:
                lcm = lcm * d // np.gcd(lcm, d)
            
            # Multiply to get integers
            coefficients = [int(f * lcm) for f in fractions]
            
            # Reduce to smallest integers
            gcd = coefficients[0]
            for c in coefficients[1:]:
                gcd = np.gcd(gcd, c)
            coefficients = [c // gcd for c in coefficients]
            
            # Create result dictionary
            result = {}
            for i, compound in enumerate(all_compounds):
                result[compound] = coefficients[i]
            
            return result
            
        except Exception as e:
            raise ValueError(f"Could not balance equation: {e}")


# ============================================================================
# MODULE 6: STOICHIOMETRY SOLVER
# ============================================================================

class StoichiometrySolver:
    """Performs mass-to-mass stoichiometry calculations"""
    
    @staticmethod
    def solve(balanced_eq: Dict[str, int], given_compound: str, 
              given_mass: float, target_compound: str) -> Dict:
        """
        Solve a mass-to-mass stoichiometry problem
        
        Args:
            balanced_eq: Dictionary of compounds to coefficients
            given_compound: Formula of the given substance
            given_mass: Mass of given substance in grams
            target_compound: Formula of target substance
            
        Returns:
            Dictionary with solution steps and final answer
        """
        # Calculate molar masses
        given_molar_mass = MolarMassCalculator.calculate(given_compound)
        target_molar_mass = MolarMassCalculator.calculate(target_compound)
        
        # Get coefficients
        given_coeff = balanced_eq[given_compound]
        target_coeff = balanced_eq[target_compound]
        
        # Step 1: Convert given mass to moles
        given_moles = given_mass / given_molar_mass
        
        # Step 2: Use mole ratio to find target moles
        target_moles = given_moles * (target_coeff / given_coeff)
        
        # Step 3: Convert target moles to mass
        target_mass = target_moles * target_molar_mass
        
        return {
            'given_molar_mass': given_molar_mass,
            'target_molar_mass': target_molar_mass,
            'given_moles': given_moles,
            'target_moles': target_moles,
            'target_mass': target_mass,
            'mole_ratio': f"{target_coeff}/{given_coeff}"
        }


# ============================================================================
# MODULE 7: MAIN INTERFACE
# ============================================================================

class StoichiometryCalculator:
    """Main interface for the stoichiometry solver"""
    
    def solve_problem(self, equation: str, given_compound: str, 
                     given_mass: float, target_compound: str, 
                     show_steps: bool = True) -> float:
        """
        Solve a complete stoichiometry problem
        
        Args:
            equation: Unbalanced chemical equation
            given_compound: Formula of given substance
            given_mass: Mass of given substance in grams
            target_compound: Formula of target substance
            show_steps: Whether to print step-by-step solution
            
        Returns:
            Mass of target substance in grams
        """
        # Parse equation
        parsed = EquationParser.parse(equation)
        reactants = parsed['reactants']
        products = parsed['products']
        
        # Balance equation
        balanced = EquationBalancer.balance(reactants, products)
        
        # Solve stoichiometry
        solution = StoichiometrySolver.solve(
            balanced, given_compound, given_mass, target_compound
        )
        
        if show_steps:
            self._print_solution(equation, balanced, given_compound, 
                               given_mass, target_compound, solution)
        
        return solution['target_mass']
    
    def _print_solution(self, equation: str, balanced: Dict[str, int],
                       given_compound: str, given_mass: float,
                       target_compound: str, solution: Dict):
        """Print formatted solution with steps"""
        print("\n" + "="*70)
        print("STOICHIOMETRY SOLUTION")
        print("="*70)
        
        # Show balanced equation
        print("\n1. BALANCED EQUATION:")
        reactants = []
        products = []
        for compound, coeff in balanced.items():
            term = f"{coeff}{compound}" if coeff > 1 else compound
            if compound in equation.split('->')[0]:
                reactants.append(term)
            else:
                products.append(term)
        print(f"   {' + '.join(reactants)} → {' + '.join(products)}")
        
        # Show given information
        print(f"\n2. GIVEN:")
        print(f"   {given_mass:.3f} g of {given_compound}")
        print(f"   Find: mass of {target_compound}")
        
        # Show molar masses
        print(f"\n3. MOLAR MASSES:")
        print(f"   {given_compound}: {solution['given_molar_mass']:.3f} g/mol")
        print(f"   {target_compound}: {solution['target_molar_mass']:.3f} g/mol")
        
        # Show calculation steps
        print(f"\n4. CALCULATION:")
        print(f"   Step 1: Convert {given_compound} to moles")
        print(f"           {given_mass:.3f} g ÷ {solution['given_molar_mass']:.3f} g/mol = {solution['given_moles']:.6f} mol")
        
        print(f"\n   Step 2: Use mole ratio ({solution['mole_ratio']})")
        print(f"           {solution['given_moles']:.6f} mol × ({solution['mole_ratio']}) = {solution['target_moles']:.6f} mol")
        
        print(f"\n   Step 3: Convert to grams of {target_compound}")
        print(f"           {solution['target_moles']:.6f} mol × {solution['target_molar_mass']:.3f} g/mol = {solution['target_mass']:.3f} g")
        
        # Final answer
        print(f"\n5. FINAL ANSWER:")
        print(f"   {solution['target_mass']:.3f} g of {target_compound}")
        print("="*70 + "\n")


# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main():
    """Main function with example problems"""
    calculator = StoichiometryCalculator()
    
    print("\nSTOICHIOMETRY SOLVER - Mass-to-Mass Problems")
    print("="*70)
    
    # Example 1
    print("\n\nEXAMPLE 1:")
    print("-" * 70)
    calculator.solve_problem(
        equation="Ca(OH)2 + HCl -> CaCl2 + H2O",
        given_compound="Ca(OH)2",
        given_mass=10.0,
        target_compound="CaCl2"
    )
    
    # Example 2
    print("\n\nEXAMPLE 2:")
    print("-" * 70)
    calculator.solve_problem(
        equation="C3H8 + O2 -> CO2 + H2O",
        given_compound="C3H8",
        given_mass=25.0,
        target_compound="CO2"
    )
    
    # Interactive mode
    print("\n" + "="*70)
    print("INTERACTIVE MODE")
    print("="*70)
    print("\nEnter your own problem (or press Enter to skip):\n")
    
    equation = input("Equation (e.g., H2 + O2 -> H2O): ").strip()
    if equation:
        given_compound = input("Given compound: ").strip()
        given_mass = float(input("Given mass (g): "))
        target_compound = input("Target compound: ").strip()
        
        calculator.solve_problem(
            equation=equation,
            given_compound=given_compound,
            given_mass=given_mass,
            target_compound=target_compound
        )


if __name__ == "__main__":
    main()