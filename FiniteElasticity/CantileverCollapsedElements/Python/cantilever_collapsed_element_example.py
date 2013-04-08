#!/usr/bin/env python

#> \file
#> \author Thiranja, Prasad, Babarenda Gamage
#> \brief This example solves a finite elasticity equation in a mesh which includes a collapsed element using openCMISS calls in python.
#>
#> \section LICENSE
#>
#> Version: MPL 1.1/GPL 2.0/LGPL 2.1
#>
#> The contents of this file are subject to the Mozilla Public License
#> Version 1.1 (the "License"); you may not use this file except in
#> compliance with the License. You may obtain a copy of the License at
#> http://www.mozilla.org/MPL/
#>
#> Software distributed under the License is distributed on an "AS IS"
#> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
#> License for the specific language governing rights and limitations
#> under the License.
#>
#> The Original Code is openCMISS
#>
#> The Initial Developer of the Original Code is University of Auckland,
#> Auckland, New Zealand and University of Oxford, Oxford, United
#> Kingdom. Portions created by the University of Auckland and University
#> of Oxford are Copyright (C) 2007 by the University of Auckland and
#> the University of Oxford. All Rights Reserved.
#>
#> Contributor(s): 
#>
#> Alternatively, the contents of this file may be used under the terms of
#> either the GNU General Public License Version 2 or later (the "GPL"), or
#> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
#> in which case the provisions of the GPL or the LGPL are applicable instead
#> of those above. if you wish to allow use of your version of this file only
#> under the terms of either the GPL or the LGPL, and not to allow others to
#> use your version of this file under the terms of the MPL, indicate your
#> decision by deleting the provisions above and replace them with the notice
#> and other provisions required by the GPL or the LGPL. if you do not delete
#> the provisions above, a recipient may use your version of this file under
#> the terms of any one of the MPL, the GPL or the LGPL.
#>

#> \example FiniteElasticity/Cantilever/src/CantileverExample.py
## Example script to solve a finite elasticity equation in a mesh which includes a collapsed element using openCMISS calls in python.
## \par Latest Builds:
## \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/FiniteElasticity/Cantilever/build-intel'>Linux Intel Build</a>
## \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/FiniteElasticity/Cantilever/build-gnu'>Linux GNU Build</a>
#<


''' About this example

This example solves a finite elasticity cantilever beam problem in which 
one of its two elements is collapsed. An additional region is created in 
order to create an OpenCMISS generated mesh whose field values will be copied 
to the collapsed mesh's geometric field. The node_mapping array specifies the 
mapping between nodes in uncollapsed generated mesh and nodes in the collapsed 
mesh.

'''

#> Main script
# Add Python bindings directory to PATH
import sys
import os

sys.path.append(os.sep.join((
    os.environ['OPENCMISS_ROOT'],'cm','bindings','python')))

# Intialise OpenCMISS
from opencmiss import CMISS

# Set problem options
collapsed = True
Hermite = True
use_pressure_basis = True
apply_load = True # For testing - solver residual should be zero if false.  
unit_scaling = False
number_of_xi = 3
number_of_global_x_elements = 2 # Need to update mapping below if changed.  
number_of_global_y_elements = 1 # Need to update mapping below if changed.  
number_of_global_z_elements = 1 # Need to update mapping below if changed.  
width = 60.0
length = 40.0
height = 40.0
density = 9.0E-4 #in g mm^-3

c1 = 1.0
number_of_load_increments = 2

# Set problem parameters
MESH_COMPONENT_1 = 1
MESH_COMPONENT_2 = 2

coordinate_system_user_number = 1
region_user_number = 1
temp_region_user_number = 2
basis_user_number = 1
collapsed_basis_user_number = 2
pressure_basis_user_number = 3
collapsed_pressure_basis_user_number = 4
generated_mesh_user_number = 1
mesh_user_number = 1
collapsed_mesh_user_number = 2
decomposition_user_number = 1
temp_decomposition_user_number = 2
geometric_field_user_number = 1
material_field_user_number = 3
dependent_field_user_number = 4
source_field_user_number = 5
equations_set_field_user_number = 6
equations_set_user_number = 1
problem_user_number = 1

if Hermite:
    interpolation_type = (
        CMISS.BasisInterpolationSpecifications.CUBIC_HERMITE)
    NumberOfGaussXi = 4
else:
    interpolation_type = (
        CMISS.BasisInterpolationSpecifications.LINEAR_LAGRANGE)
    NumberOfGaussXi = 2

number_of_elements = (number_of_global_x_elements *
                      number_of_global_y_elements *
                      number_of_global_z_elements)

if collapsed:
    number_of_nodes = 10
    node_mapping = [1, 2, 3, 4, 5, 6, 7, 8, 10, 11]
    boundary_condition_nodes = [1, 4, 7, 9] # Fixed end of cantilever.  
else:
    number_of_nodes = 12
    node_mapping = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    boundary_condition_nodes = [1, 4, 7, 10] # Fixed end of cantilever.  

if unit_scaling:
    scaling_type = CMISS.FieldScalingTypes.UNIT
else:
    scaling_type = CMISS.FieldScalingTypes.ARITHMETIC_MEAN

if apply_load:
    gravity=[0.0, 0.0, -9.81] #in m s^-2
else:
    gravity=[0.0, 0.0, 0.0] #in m s^-2

# Set all diganostic levels on for testing

# CMISS.DiagnosticsSetOn(CMISS.DiagnosticTypes.All,
#                        [1, 2, 3, 4, 5],
#                        'Diagnostics',
#                        ['DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE'])

# Get the number of computational nodes and this computational node number.  
numberOfComputationalNodes = CMISS.ComputationalNumberOfNodesGet()
computationalNodeNumber = CMISS.ComputationalNodeNumberGet()

# Create a 3D rectangular cartesian coordinate system.  
coordinateSystem = CMISS.CoordinateSystem()
coordinateSystem.CreateStart(coordinate_system_user_number)
coordinateSystem.DimensionSet(3)
coordinateSystem.CreateFinish()

# Create a region and assign the coordinate system to the region.  
region = CMISS.Region()
region.CreateStart(region_user_number, CMISS.WorldRegion)
region.LabelSet('Region')
region.CoordinateSystemSet(coordinateSystem)
region.CreateFinish()

# Temporary region for creating a regular rectangular generated mesh whoses 
# geometric field will be used to populate the collapsed mesh's geometric 
# field.  
temp_region = CMISS.Region()
temp_region.CreateStart(temp_region_user_number, CMISS.WorldRegion)
temp_region.LabelSet('temp_Region')
temp_region.CoordinateSystemSet(coordinateSystem)
temp_region.CreateFinish()

# Define basis
basis = CMISS.Basis()
basis.CreateStart(basis_user_number)
basis.TypeSet(CMISS.BasisTypes.LAGRANGE_HERMITE_TP)
basis.NumberOfXiSet(number_of_xi)
basis.InterpolationXiSet([interpolation_type]*number_of_xi)
basis.QuadratureNumberOfGaussXiSet([NumberOfGaussXi]*number_of_xi)
basis.CreateFinish()

if collapsed:
    # Define collapsed basis and specify which xi directions and collapsed 
    # and which are not using CollapsedXiSet.  
    # XI_COLLAPSED = 1      # The Xi direction is collapsed
    # COLLAPSED_AT_XI0 = 2  # The Xi direction at the xi=0 end of this Xi 
    #                         direction is collapsed
    # COLLAPSED_AT_XI1 = 3  # The Xi direction at the xi=1 end of this Xi 
    #                         direction is collapsed
    # NOT_COLLAPSED = 4     # The Xi direction is not collapsed
    collapsed_basis = CMISS.Basis()
    collapsed_basis.CreateStart(collapsed_basis_user_number)
    collapsed_basis.TypeSet(CMISS.BasisTypes.LAGRANGE_HERMITE_TP)
    collapsed_basis.NumberOfXiSet(number_of_xi)
    collapsed_basis.InterpolationXiSet([interpolation_type]*number_of_xi)
    collapsed_basis.CollapsedXiSet([CMISS.BasisXiCollapse.COLLAPSED_AT_XI1,
                                    CMISS.BasisXiCollapse.NOT_COLLAPSED,
                                    CMISS.BasisXiCollapse.XI_COLLAPSED])
    collapsed_basis.QuadratureNumberOfGaussXiSet(
        [NumberOfGaussXi]*number_of_xi)
    collapsed_basis.CreateFinish()

if(use_pressure_basis):
    # Define pressure basis
    pressure_basis = CMISS.Basis()
    pressure_basis.CreateStart(pressure_basis_user_number)
    pressure_basis.TypeSet(CMISS.BasisTypes.LAGRANGE_HERMITE_TP)
    pressure_basis.NumberOfXiSet(number_of_xi)
    pressure_basis.InterpolationXiSet(
        [CMISS.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*number_of_xi)
    pressure_basis.QuadratureNumberOfGaussXiSet(
        [NumberOfGaussXi]*number_of_xi)
    pressure_basis.CreateFinish()

    if collapsed:
        collapsed_pressure_basis = CMISS.Basis()
        collapsed_pressure_basis.CreateStart(
            collapsed_pressure_basis_user_number)
        collapsed_pressure_basis.TypeSet(CMISS.BasisTypes.LAGRANGE_HERMITE_TP)
        collapsed_pressure_basis.NumberOfXiSet(number_of_xi)
        collapsed_pressure_basis.InterpolationXiSet(
            [CMISS.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*
             number_of_xi)
        collapsed_pressure_basis.CollapsedXiSet(
            [CMISS.BasisXiCollapse.COLLAPSED_AT_XI1,
             CMISS.BasisXiCollapse.NOT_COLLAPSED,
             CMISS.BasisXiCollapse.XI_COLLAPSED])
        collapsed_pressure_basis.QuadratureNumberOfGaussXiSet(
            [NumberOfGaussXi]*number_of_xi)
        collapsed_pressure_basis.CreateFinish()

mesh = CMISS.Mesh() # Create the OpenCMISS mesh object.  
mesh.CreateStart(collapsed_mesh_user_number, region, number_of_xi)
mesh.NumberOfComponentsSet(1 + use_pressure_basis)
mesh.NumberOfElementsSet(number_of_elements)

nodes = CMISS.Nodes()
nodes.CreateStart(region, number_of_nodes)
nodes.CreateFinish()

elements = CMISS.MeshElements()
elements.CreateStart(mesh, MESH_COMPONENT_1, basis)
if collapsed:
    elements.NodesSet(1, [1, 2, 4, 5, 7, 8, 9, 10])
    elements.BasisSet(2, collapsed_basis)
    elements.NodesSet(2, [2, 3, 5, 6, 8, 10])
    collapsed_nodes = [3, 6]
else:
    elements.NodesSet(1, [1, 2, 4, 5, 7, 8, 10, 11])
    elements.NodesSet(2, [2, 3, 5, 6, 8, 9, 11, 12])
    collapsed_nodes = []
elements.CreateFinish()

if(use_pressure_basis):
    pressure_elements = CMISS.MeshElements()
    pressure_elements.CreateStart(mesh, MESH_COMPONENT_2, pressure_basis)
    if collapsed:
        pressure_elements.NodesSet(1, [1, 2, 4, 5, 7, 8, 9, 10])
        pressure_elements.BasisSet(2, collapsed_pressure_basis)
        pressure_elements.NodesSet(2, [2, 3, 5, 6, 8, 10])
        collapsed_nodes = [3, 6]
    else:
        pressure_elements.NodesSet(1, [1, 2, 4, 5, 7, 8, 10, 11])
        pressure_elements.NodesSet(2, [2, 3, 5, 6, 8, 9, 11, 12])
        collapsed_nodes = []
    pressure_elements.CreateFinish()

mesh.CreateFinish()

# Start the creation of a generated mesh in the region
generatedMesh = CMISS.GeneratedMesh()
generatedMesh.CreateStart(generated_mesh_user_number, temp_region)
generatedMesh.TypeSet(CMISS.GeneratedMeshTypes.REGULAR)
if(use_pressure_basis):
    generatedMesh.BasisSet([basis, pressure_basis])
else:
    generatedMesh.BasisSet([basis])
generatedMesh.ExtentSet([width, length, height])
generatedMesh.NumberOfElementsSet([number_of_global_x_elements,
                                   number_of_global_y_elements,
                                   number_of_global_z_elements])
# Finish the creation of a generated mesh in the region
temp_mesh = CMISS.Mesh()
generatedMesh.CreateFinish(mesh_user_number, temp_mesh)

# Create a decomposition for the mesh
decomposition = CMISS.Decomposition()
decomposition.CreateStart(decomposition_user_number, mesh)
decomposition.TypeSet(CMISS.DecompositionTypes.CALCULATED)
decomposition.NumberOfDomainsSet(numberOfComputationalNodes)
decomposition.CalculateFacesSet(True)
decomposition.CreateFinish()

# Create a decomposition for the mesh
temp_decomposition = CMISS.Decomposition()
temp_decomposition.CreateStart(temp_decomposition_user_number, temp_mesh)
temp_decomposition.TypeSet(CMISS.DecompositionTypes.CALCULATED)
temp_decomposition.NumberOfDomainsSet(numberOfComputationalNodes)
temp_decomposition.CalculateFacesSet(True)
temp_decomposition.CreateFinish()

# Create a field for the geometry
temp_geometric_field = CMISS.Field()
temp_geometric_field.CreateStart(geometric_field_user_number, temp_region)
temp_geometric_field.MeshDecompositionSet(temp_decomposition)
temp_geometric_field.TypeSet(CMISS.FieldTypes.GEOMETRIC)
temp_geometric_field.VariableLabelSet(
    CMISS.FieldVariableTypes.U,'temp_Geometry')
temp_geometric_field.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U,1,1)
temp_geometric_field.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U,2,1)
temp_geometric_field.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U,3,1)
temp_geometric_field.ScalingTypeSet(scaling_type)
temp_geometric_field.CreateFinish()

# Update the geometric field parameters from generated mesh
generatedMesh.GeometricParametersCalculate(temp_geometric_field)

# Create a field for the geometry
geometric_field = CMISS.Field()
geometric_field.CreateStart(geometric_field_user_number, region)
geometric_field.LabelSet('Geometry')
geometric_field.MeshDecompositionSet(decomposition)
geometric_field.TypeSet(CMISS.FieldTypes.GEOMETRIC)
geometric_field.VariableLabelSet(CMISS.FieldVariableTypes.U,'Geometry')
geometric_field.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U,1,1)
geometric_field.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U,2,1)
geometric_field.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U,3,1)
geometric_field.ScalingTypeSet(scaling_type)
geometric_field.CreateFinish()

# Update the geometric field parameters from generated mesh.  
for collapsed_mesh_node in range(1, number_of_nodes + 1):
    uncollapsed_mesh_node = node_mapping[collapsed_mesh_node - 1]
    # Copy nodal geometric values.  
    for component in [1, 2, 3]:
        value = temp_geometric_field.ParameterSetGetNodeDP(
            CMISS.FieldVariableTypes.U, 
            CMISS.FieldParameterSetTypes.VALUES,
            1, 1, uncollapsed_mesh_node, component)
        geometric_field.ParameterSetUpdateNodeDP(
            CMISS.FieldVariableTypes.U,
            CMISS.FieldParameterSetTypes.VALUES,
            1, 1, collapsed_mesh_node, component, value)
    # Copy nodal derivative values.  
    if collapsed_mesh_node in collapsed_nodes:
      component = 2
      collapsed_derivative = 2
      uncollapsed_derivative = 3
      value = temp_geometric_field.ParameterSetGetNodeDP(
          CMISS.FieldVariableTypes.U,
          CMISS.FieldParameterSetTypes.VALUES,
          1, uncollapsed_derivative, uncollapsed_mesh_node, component)
      geometric_field.ParameterSetUpdateNodeDP(
          CMISS.FieldVariableTypes.U,
          CMISS.FieldParameterSetTypes.VALUES,
          1, collapsed_derivative, collapsed_mesh_node, component, value)
    else:
        for component in [1, 2, 3]:
            if Hermite:
                for derivative in range(2, 9):
                    value = temp_geometric_field.ParameterSetGetNodeDP(
                        CMISS.FieldVariableTypes.U,
                        CMISS.FieldParameterSetTypes.VALUES,
                        1, derivative, uncollapsed_mesh_node, component)
                    geometric_field.ParameterSetUpdateNodeDP(
                        CMISS.FieldVariableTypes.U,
                        CMISS.FieldParameterSetTypes.VALUES,
                        1, derivative, collapsed_mesh_node, component, value)

geometric_field.ParameterSetUpdateStart(
    CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES)
geometric_field.ParameterSetUpdateFinish(
    CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES)

# Create the equations_set
equationsSetField = CMISS.Field()
equationsSet = CMISS.EquationsSet()
equationsSet.CreateStart(equations_set_user_number, region, geometric_field,
                         CMISS.EquationsSetClasses.ELASTICITY,
                         CMISS.EquationsSetTypes.FINITE_ELASTICITY,
                         CMISS.EquationsSetSubtypes.MOONEY_RIVLIN,
                         equations_set_field_user_number, equationsSetField)
equationsSet.CreateFinish()

# Create the dependent field
dependentField = CMISS.Field()
equationsSet.DependentCreateStart(dependent_field_user_number, dependentField)
dependentField.LabelSet('Dependent')
dependentField.VariableLabelSet(CMISS.FieldVariableTypes.U, 'Dependent')
if(use_pressure_basis):
    # Use pressure mesh component 
    dependentField.ComponentMeshComponentSet(
        CMISS.FieldVariableTypes.U, 4, 2)
    dependentField.ComponentMeshComponentSet(
        CMISS.FieldVariableTypes.DELUDELN, 4, 2)
    # Set the pressure to be nodally based and use the second mesh component.  
    if Hermite:
        dependentField.ComponentInterpolationSet(
            CMISS.FieldVariableTypes.U, 4,
            CMISS.FieldInterpolationTypes.NODE_BASED)
        dependentField.ComponentInterpolationSet(
            CMISS.FieldVariableTypes.DELUDELN, 4,
            CMISS.FieldInterpolationTypes.NODE_BASED)
else:
    dependentField.ComponentInterpolationSet(
        CMISS.FieldVariableTypes.U, 4,
        CMISS.FieldInterpolationTypes.ELEMENT_BASED)
    dependentField.ComponentInterpolationSet(
        CMISS.FieldVariableTypes.DELUDELN, 4,
        CMISS.FieldInterpolationTypes.ELEMENT_BASED)
dependentField.ScalingTypeSet(scaling_type)
equationsSet.DependentCreateFinish()


# Initialise dependent field parmeters from the geometric field parameters.  
for component in [1, 2, 3]:
    CMISS.Field.ParametersToFieldParametersComponentCopy(
        geometric_field, CMISS.FieldVariableTypes.U,
        CMISS.FieldParameterSetTypes.VALUES, component,
        dependentField, CMISS.FieldVariableTypes.U,
        CMISS.FieldParameterSetTypes.VALUES, component)
CMISS.Field.ComponentValuesInitialiseDP(
    dependentField,CMISS.FieldVariableTypes.U,
    CMISS.FieldParameterSetTypes.VALUES, 4, -c1)

material_field = CMISS.Field()
material_field.CreateStart(material_field_user_number, region)
material_field.LabelSet('Material')
material_field.MeshDecompositionSet(decomposition)
material_field.TypeSet(CMISS.FieldTypes.MATERIAL)
material_field.GeometricFieldSet(geometric_field)
material_field.NumberOfVariablesSet(2)
material_field.VariableTypesSet([CMISS.FieldVariableTypes.U,
                                  CMISS.FieldVariableTypes.V])
material_field.NumberOfComponentsSet(CMISS.FieldVariableTypes.U, 2)
material_field.NumberOfComponentsSet(CMISS.FieldVariableTypes.V, 1)
material_field.ComponentInterpolationSet(CMISS.FieldVariableTypes.U, 1,
    CMISS.FieldInterpolationTypes.ELEMENT_BASED)
material_field.ComponentInterpolationSet(CMISS.FieldVariableTypes.U, 2,
    CMISS.FieldInterpolationTypes.ELEMENT_BASED)
material_field.ComponentInterpolationSet(CMISS.FieldVariableTypes.V, 1,
    CMISS.FieldInterpolationTypes.ELEMENT_BASED)
material_field.VariableLabelSet(CMISS.FieldVariableTypes.U, 'Material')
material_field.VariableLabelSet(CMISS.FieldVariableTypes.V, 'Density')
material_field.ScalingTypeSet(scaling_type)
material_field.CreateFinish()

# Create the material field
equationsSet.MaterialsCreateStart(material_field_user_number, material_field)
equationsSet.MaterialsCreateFinish()

# Set Mooney-Rivlin constants c10 and c01 respectively.
material_field.ComponentValuesInitialiseDP(
    CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES,
    1 , c1)
material_field.ComponentValuesInitialiseDP(
    CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES,
    2, 0.0)
material_field.ComponentValuesInitialiseDP(
    CMISS.FieldVariableTypes.V, CMISS.FieldParameterSetTypes.VALUES,
    1, density)

source_field = CMISS.Field()
source_field.CreateStart(source_field_user_number, region)
source_field.LabelSet('Gravity')
source_field.MeshDecompositionSet(decomposition)
source_field.TypeSet(CMISS.FieldTypes.GENERAL)
source_field.GeometricFieldSet(geometric_field)
for component in [1, 2, 3]:
    source_field.ComponentInterpolationSet(
        CMISS.FieldVariableTypes.U, component,
        CMISS.FieldInterpolationTypes.ELEMENT_BASED)
source_field.ScalingTypeSet(scaling_type)
source_field.CreateFinish()

#Create the source field with the gravity vector
equationsSet.SourceCreateStart(source_field_user_number, source_field)
equationsSet.SourceCreateFinish()

#Set the gravity vector component values
for component in [1, 2, 3]:
    source_field.ComponentValuesInitialiseDP(
        CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES,
        component, gravity[component - 1])

# Create equations
equations = CMISS.Equations()
equationsSet.EquationsCreateStart(equations)
equations.sparsityType = CMISS.EquationsSparsityTypes.SPARSE
equations.outputType = CMISS.EquationsOutputTypes.NONE
equationsSet.EquationsCreateFinish()

# Define the problem
problem = CMISS.Problem()
problem.CreateStart(problem_user_number)
problem.SpecificationSet(CMISS.ProblemClasses.ELASTICITY,
                         CMISS.ProblemTypes.FINITE_ELASTICITY,
                         CMISS.ProblemSubTypes.NONE)
problem.CreateFinish()

# Create the problem control loop
problem.ControlLoopCreateStart()
controlLoop = CMISS.ControlLoop()
problem.ControlLoopGet([CMISS.ControlLoopIdentifiers.NODE], controlLoop)
controlLoop.MaximumIterationsSet(number_of_load_increments)
problem.ControlLoopCreateFinish()

# Create problem solver
nonLinearSolver = CMISS.Solver()
linearSolver = CMISS.Solver()
problem.SolversCreateStart()
problem.SolverGet([CMISS.ControlLoopIdentifiers.NODE], 1, nonLinearSolver)
nonLinearSolver.OutputTypeSet(CMISS.SolverOutputTypes.PROGRESS)
nonLinearSolver.NewtonJacobianCalculationTypeSet(
    CMISS.JacobianCalculationTypes.EQUATIONS)
nonLinearSolver.NewtonLinearSolverGet(linearSolver)
linearSolver.LinearTypeSet(CMISS.LinearSolverTypes.DIRECT)
#linearSolver.LibraryTypeSet(CMISS.SolverLibraries.LAPACK)
problem.SolversCreateFinish()

# Create solver equations and add equations set to solver equations
solver = CMISS.Solver()
solverEquations = CMISS.SolverEquations()
problem.SolverEquationsCreateStart()
problem.SolverGet([CMISS.ControlLoopIdentifiers.NODE], 1, solver)
solver.SolverEquationsGet(solverEquations)
solverEquations.SparsityTypeSet(CMISS.SolverEquationsSparsityTypes.SPARSE)
equationsSetIndex = solverEquations.EquationsSetAdd(equationsSet)
problem.SolverEquationsCreateFinish()

# Prescribe boundary conditions (absolute nodal parameters)
boundaryConditions = CMISS.BoundaryConditions()
solverEquations.BoundaryConditionsCreateStart(boundaryConditions)

for node in boundary_condition_nodes:
    for component in [1, 2, 3]:
        boundaryConditions.AddNode(
            dependentField, CMISS.FieldVariableTypes.U, 1, 1, 
            node, component, CMISS.BoundaryConditionsTypes.FIXED, 0.0)
        if Hermite:
            for derivative in range(3,9):
                boundaryConditions.AddNode(
                    dependentField,CMISS.FieldVariableTypes.U,1 , derivative,
                    node,component,CMISS.BoundaryConditionsTypes.FIXED,0.0)
solverEquations.BoundaryConditionsCreateFinish()

# Solve the problem
problem.Solve()

# Export results
fields = CMISS.Fields()
fields.CreateRegion(region)
fields.NodesExport('Cantilever', 'FORTRAN')
fields.ElementsExport('Cantilever', 'FORTRAN')
fields.Finalise()

print 'Program successfully completed.'
