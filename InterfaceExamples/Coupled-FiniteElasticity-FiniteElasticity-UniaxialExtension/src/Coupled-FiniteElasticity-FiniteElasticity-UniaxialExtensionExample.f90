!> \file
!> \author Chris Bradley
!> \brief This is an example program which solves a weakly coupled finite-elasticity-finite-elasticity problem in two regions using openCMISS calls.
!>
!> \section LICENSE
!>
!> Version: MPL 1.1/GPL 2.0/LGPL 2.1
!>
!> The contents of this file are subject to the Mozilla Public License
!> Version 1.1 (the "License"); you may not use this file except in
!> compliance with the License. You may obtain a copy of the License at
!> http://www.mozilla.org/MPL/
!>
!> Software distributed under the License is distributed on an "AS IS"
!> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
!> License for the specific language governing rights and limitations
!> under the License.
!>
!> The Original Code is OpenCMISS
!>
!> The Initial Developer of the Original Code is University of Auckland,
!> Auckland, New Zealand and University of Oxford, Oxford, United
!> Kingdom. Portions created by the University of Auckland and University
!> of Oxford are Copyright (C) 2007 by the University of Auckland and
!> the University of Oxford. All Rights Reserved.
!>
!> Contributor(s):Lukas Moj, Thiranja Prasad Babarenda Gamage
!>
!> Alternatively, the contents of this file may be used under the terms of
!> either the GNU General Public License Version 2 or later (the "GPL"), or
!> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
!> in which case the provisions of the GPL or the LGPL are applicable instead
!> of those above. If you wish to allow use of your version of this file only
!> under the terms of either the GPL or the LGPL, and not to allow others to
!> use your version of this file under the terms of the MPL, indicate your
!> decision by deleting the provisions above and replace them with the notice
!> and other provisions required by the GPL or the LGPL. If you do not delete
!> the provisions above, a recipient may use your version of this file under
!> the terms of any one of the MPL, the GPL or the LGPL.
!>

!> \example InterfaceExamples/Coupled-FiniteElasticity-FiniteElasticity-PenaltyMethod/src/Coupled-FiniteElasticity-FiniteElasticity-UniaxialExtensionExample.f90
!! Example program which sets up a field in two regions using OpenCMISS calls.
!! \par Latest Builds:
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/InterfaceExamples/Coupled-FiniteElasticity-FiniteElasticity/build-intel'>Linux Intel Build</a>
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/InterfaceExamples/Coupled-FiniteElasticity-FiniteElasticity/build-gnu'>Linux GNU Build</a>
!<

!> Main program
PROGRAM COUPLEDFINITEELASTICITYFINITEELASTICITY

  USE OPENCMISS
  USE MPI
  
#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

  !Test program parameters
 
  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystem1UserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystem2UserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: Region1UserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: Region2UserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: DisplacementBasis1UserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: DisplacementBasis2UserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: HydrostaticPressureBasis1UserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: HydrostaticPressureBasis2UserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: InterfaceBasisUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: InterfaceMappingBasisUserNumber=6
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMesh1UserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMesh2UserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: InterfaceGeneratedMeshUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: Mesh1UserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: Mesh2UserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: InterfaceMeshUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: Decomposition1UserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: Decomposition2UserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: InterfaceDecompositionUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: GeometricField1UserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: GeometricField2UserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: InterfaceGeometricFieldUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: FibreField1UserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: FibreField2UserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: MaterialField1UserNumber=6
  INTEGER(CMISSIntg), PARAMETER :: MaterialField2UserNumber=7
  INTEGER(CMISSIntg), PARAMETER :: DependentField1UserNumber=8
  INTEGER(CMISSIntg), PARAMETER :: DependentField2UserNumber=9
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetField1UserNumber=10
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetField2UserNumber=11
  INTEGER(CMISSIntg), PARAMETER :: LagrangeFieldUserNumber=12
  INTEGER(CMISSIntg), PARAMETER :: PenaltyFieldUserNumber=13
  INTEGER(CMISSIntg), PARAMETER :: GeometricField1NumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: GeometricField1NumberOfComponents=3
  INTEGER(CMISSIntg), PARAMETER :: GeometricField2NumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: GeometricField2NumberOfComponents=3
  INTEGER(CMISSIntg), PARAMETER :: FibreField1NumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: FibreField1NumberOfComponents=3
  INTEGER(CMISSIntg), PARAMETER :: FibreField2NumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: FibreField2NumberOfComponents=3
  INTEGER(CMISSIntg), PARAMETER :: MaterialField1NumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: MaterialField2NumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: DependentField1NumberOfVariables=2
  INTEGER(CMISSIntg), PARAMETER :: DependentField2NumberOfVariables=2
  INTEGER(CMISSIntg), PARAMETER :: InterfaceUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: InterfaceConditionUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: EquationsSet1UserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: EquationsSet2UserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: CoupledProblemUserNumber=1

  !Program types
  
  !Program variables

  INTEGER(CMISSIntg) :: NumberGlobalXElements,NumberGlobalYElements,NumberGlobalZElements
  INTEGER(CMISSIntg) :: DisplacementInterpolationType,PressureInterpolationType,NumberOfGaussXi,NumberOfNodeXi

  INTEGER(CMISSIntg) :: EquationsSet1Index,EquationsSet2Index
  INTEGER(CMISSIntg) :: InterfaceConditionIndex
  INTEGER(CMISSIntg) :: Mesh1Index,Mesh2Index
  INTEGER(CMISSIntg) :: MaterialField1NumberOfComponents,MaterialField2NumberOfComponents
  INTEGER(CMISSIntg) :: DependentField1NumberOfComponents,DependentField2NumberOfComponents
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,ComputationalNodeNumber
  INTEGER(CMISSIntg) :: y_element_idx,z_element_idx,mesh_local_y_node,mesh_local_z_node
  INTEGER(CMISSIntg) :: FrontNodeGroup(4),BackNodeGroup(4),Side4NodeGroup(4),Side6NodeGroup(6)
  INTEGER(CMISSIntg) :: Bottom4NodeGroup(4),Bottom6NodeGroup(6),node_idx,node,component,LoadSteps
  REAL(CMISSDP) :: XI3(3), Height,Width,Length
  LOGICAL :: HydrostaticPressureBasis,Uncoupled,ForceBoundaryConditions,Incompressible1,Incompressible2,PenaltyMethod

  !CMISS variables

  TYPE(CMISSBasisType) :: DisplacementBasis1,DisplacementBasis2,HydrostaticPressureBasis1,HydrostaticPressureBasis2, &
    & InterfaceBasis,InterfaceMappingBasis
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditions
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem1,CoordinateSystem2,WorldCoordinateSystem
  TYPE(CMISSDecompositionType) :: Decomposition1,Decomposition2,InterfaceDecomposition
  TYPE(CMISSEquationsType) :: Equations1,Equations2
  TYPE(CMISSEquationsSetType) :: EquationsSet1,EquationsSet2
  TYPE(CMISSFieldType) :: GeometricField1,GeometricField2,InterfaceGeometricField,DependentField1,DependentField2,LagrangeField, &
    & EquationsSetField1,EquationsSetField2,PenaltyField
  TYPE(CMISSFieldType) :: FibreField1,FibreField2,MaterialField1,MaterialField2
  TYPE(CMISSFieldsType) :: Fields1,Fields2,Fields3
  TYPE(CMISSGeneratedMeshType) :: GeneratedMesh1,GeneratedMesh2,InterfaceGeneratedMesh
  TYPE(CMISSMeshType) :: Mesh1,Mesh2,InterfaceMesh
  TYPE(CMISSInterfaceMeshConnectivityType) :: InterfaceMeshConnectivity
  TYPE(CMISSInterfaceType) :: Interface
  TYPE(CMISSInterfaceConditionType) :: InterfaceCondition
  TYPE(CMISSInterfaceEquationsType) :: InterfaceEquations
  TYPE(CMISSProblemType) :: CoupledProblem
  TYPE(CMISSRegionType) :: Region1,Region2,WorldRegion
  TYPE(CMISSSolverType) :: NonLinearSolver,LinearSolver
  TYPE(CMISSSolverEquationsType) :: NonLinearSolverEquations
  TYPE(CMISSControlLoopType) :: ControlLoop
  
#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif
  
  !Generic CMISS variables
  
  INTEGER(CMISSIntg) :: Err
  
#ifdef WIN32
  !Initialise QuickWin
  QUICKWIN_WINDOW_CONFIG%TITLE="General Output" !Window title
  QUICKWIN_WINDOW_CONFIG%NUMTEXTROWS=-1 !Max possible number of rows
  QUICKWIN_WINDOW_CONFIG%MODE=QWIN$SCROLLDOWN
  !Set the window parameters
  QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
  !If attempt fails set with system estimated values
  IF(.NOT.QUICKWIN_STATUS) QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
#endif

  PenaltyMethod=.FALSE.
  ForceBoundaryConditions=.TRUE. 
  Incompressible1=.FALSE.
  Incompressible2=.FALSE. 
  HydrostaticPressureBasis=.FALSE.
  Uncoupled=.FALSE. 
  IF(Uncoupled.EQV..TRUE.) THEN
    NumberGlobalXElements=2
    Width=2.0_CMISSDP
  ELSE
    NumberGlobalXElements=1
    Width=1.0_CMISSDP
  ENDIF
  Height=1.0_CMISSDP
  Length=1.0_CMISSDP
  NumberGlobalYElements=1
  NumberGlobalZElements=1
  !DisplacementInterpolationType = CMISS_BASIS_LINEAR_LAGRANGE_INTERPOLATION
  DisplacementInterpolationType = CMISS_BASIS_CUBIC_HERMITE_INTERPOLATION
  PressureInterpolationType=CMISS_BASIS_LINEAR_LAGRANGE_INTERPOLATION
  SELECT CASE(DisplacementInterpolationType)
  CASE(CMISS_BASIS_LINEAR_LAGRANGE_INTERPOLATION)
    NumberOfGaussXi=2
  CASE(CMISS_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION)
    NumberOfGaussXi=3
  CASE(CMISS_BASIS_CUBIC_LAGRANGE_INTERPOLATION,CMISS_BASIS_CUBIC_HERMITE_INTERPOLATION)
    NumberOfGaussXi=4
  CASE DEFAULT
    CALL HANDLE_ERROR("Invalid interpolation type.")
  END SELECT
  IF (Incompressible1.eqv..TRUE.) THEN
    MaterialField1NumberOfComponents=2
    DependentField1NumberOfComponents=4
  ELSE
    MaterialField1NumberOfComponents=3
    DependentField1NumberOfComponents=3
  ENDIF

  IF (Incompressible2.eqv..TRUE.) THEN
    MaterialField2NumberOfComponents=2
    DependentField2NumberOfComponents=4
  ELSE
    MaterialField2NumberOfComponents=3
    DependentField2NumberOfComponents=3
  ENDIF
  LoadSteps=1

  !Intialise OpenCMISS
  CALL CMISSInitialise(WorldCoordinateSystem,WorldRegion,Err)

  !Set error handling mode
  CALL CMISSErrorHandlingModeSet(CMISS_ERRORS_TRAP_ERROR,Err)
 
  !Set diganostics for testing
  !CALL CMISSDiagnosticsSetOn(CMISS_ALL_DIAG_TYPE,[1,2,3,4,5],"Diagnostics",["SOLVER_MAPPING_CALCULATE         ", &
  !  & "SOLVER_MATRIX_STRUCTURE_CALCULATE"],Err)
  
  !Get the computational nodes information
  CALL CMISSComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL CMISSComputationalNodeNumberGet(ComputationalNodeNumber,Err)
  
  !Start the creation of a new RC coordinate system for the first region
  PRINT *, ' == >> CREATING COORDINATE SYSTEM(1) << == '
  CALL CMISSCoordinateSystem_Initialise(CoordinateSystem1,Err)
  CALL CMISSCoordinateSystem_CreateStart(CoordinateSystem1UserNumber,CoordinateSystem1,Err)
  !Set the coordinate system to be 3D
  CALL CMISSCoordinateSystem_DimensionSet(CoordinateSystem1,3,Err)
  CALL CMISSCoordinateSystem_OriginSet(CoordinateSystem1,(/0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP/),Err)
  !Finish the creation of the coordinate system
  CALL CMISSCoordinateSystem_CreateFinish(CoordinateSystem1,Err)   

  !Start the creation of a new RC coordinate system for the second region
  PRINT *, ' == >> CREATING COORDINATE SYSTEM(2) << == '
  CALL CMISSCoordinateSystem_Initialise(CoordinateSystem2,Err)
  CALL CMISSCoordinateSystem_CreateStart(CoordinateSystem2UserNumber,CoordinateSystem2,Err)
  !Set the coordinate system to be 3D
  CALL CMISSCoordinateSystem_DimensionSet(CoordinateSystem2,3,Err)
  CALL CMISSCoordinateSystem_OriginSet(CoordinateSystem2,(/0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP/),Err)
  !Finish the creation of the coordinate system
  CALL CMISSCoordinateSystem_CreateFinish(CoordinateSystem2,Err) 
  
  !Start the creation of the first region
  PRINT *, ' == >> CREATING REGION(1) << == '
  CALL CMISSRegion_Initialise(Region1,Err)
  CALL CMISSRegion_CreateStart(Region1UserNumber,WorldRegion,Region1,Err)
  CALL CMISSRegion_LabelSet(Region1,"Region1",Err)
  !Set the regions coordinate system to the 1st RC coordinate system that we have created
  CALL CMISSRegion_CoordinateSystemSet(Region1,CoordinateSystem1,Err)
  !Finish the creation of the first region
  CALL CMISSRegion_CreateFinish(Region1,Err)

  !Start the creation of the second region
  PRINT *, ' == >> CREATING REGION(2) << == '
  CALL CMISSRegion_Initialise(Region2,Err)
  CALL CMISSRegion_CreateStart(Region2UserNumber,WorldRegion,Region2,Err)
  CALL CMISSRegion_LabelSet(Region2,"Region2",Err)
  !Set the regions coordinate system to the 2nd RC coordinate system that we have created
  CALL CMISSRegion_CoordinateSystemSet(Region2,CoordinateSystem2,Err)
  !Finish the creation of the second region
  CALL CMISSRegion_CreateFinish(Region2,Err) 

  !Start the creation of a basis to describe displacement
  PRINT *, ' == >> CREATING DISPLACEMENT BASIS(1) << == '
  CALL CMISSBasis_Initialise(DisplacementBasis1,Err)
  CALL CMISSBasis_CreateStart(DisplacementBasis1UserNumber,DisplacementBasis1,Err)
  CALL CMISSBasis_NumberOfXiSet(DisplacementBasis1,3,Err)
  CALL CMISSBasis_InterpolationXiSet(DisplacementBasis1,[DisplacementInterpolationType,DisplacementInterpolationType, &
    & DisplacementInterpolationType],Err)
  CALL CMISSBasis_QuadratureNumberOfGaussXiSet(DisplacementBasis1,[NumberOfGaussXi,NumberOfGaussXi, &
    & NumberOfGaussXi],Err)
  !Finish the creation of the basis
  CALL CMISSBasis_CreateFinish(DisplacementBasis1,Err) 

  !Start the creation of a basis to describe hydrostatic pressure
  IF (HydrostaticPressureBasis.eqv..TRUE.) THEN
    PRINT *, ' == >> CREATING PRESSURE BASIS(1) << == '
    CALL CMISSBasis_Initialise(HydrostaticPressureBasis1,Err)
    CALL CMISSBasis_CreateStart(HydrostaticPressureBasis1UserNumber,HydrostaticPressureBasis1,Err)
    CALL CMISSBasis_NumberOfXiSet(HydrostaticPressureBasis1,3,Err)
    CALL CMISSBasis_InterpolationXiSet(HydrostaticPressureBasis1,[PressureInterpolationType,PressureInterpolationType, &
      & PressureInterpolationType],Err)
    CALL CMISSBasis_QuadratureNumberOfGaussXiSet(HydrostaticPressureBasis1,[NumberOfGaussXi,NumberOfGaussXi, &
      & NumberOfGaussXi],Err)
    !Finish the creation of the pressure basis
    CALL CMISSBasis_CreateFinish(HydrostaticPressureBasis1,Err)
  ENDIF
   
  !Start the creation of a basis to describe displacement
  PRINT *, ' == >> CREATING DISPLACEMENT BASIS(2) << == '
  CALL CMISSBasis_Initialise(DisplacementBasis2,Err)
  CALL CMISSBasis_CreateStart(DisplacementBasis2UserNumber,DisplacementBasis2,Err)
  CALL CMISSBasis_NumberOfXiSet(DisplacementBasis2,3,Err)
  CALL CMISSBasis_InterpolationXiSet(DisplacementBasis2,[DisplacementInterpolationType,DisplacementInterpolationType, &
    & DisplacementInterpolationType],Err)
  CALL CMISSBasis_QuadratureNumberOfGaussXiSet(DisplacementBasis2,[NumberOfGaussXi,NumberOfGaussXi, &
    & NumberOfGaussXi],Err)
  !Finish the creation of the basis
  CALL CMISSBasis_CreateFinish(DisplacementBasis2,Err)

  IF (HydrostaticPressureBasis.eqv..TRUE.) THEN
    PRINT *, ' == >> CREATING PRESSURE BASIS(2) << == '
    CALL CMISSBasis_Initialise(HydrostaticPressureBasis2,Err)
    CALL CMISSBasis_CreateStart(HydrostaticPressureBasis2UserNumber,HydrostaticPressureBasis2,Err)
    CALL CMISSBasis_NumberOfXiSet(HydrostaticPressureBasis2,3,Err)
    CALL CMISSBasis_InterpolationXiSet(HydrostaticPressureBasis2,[PressureInterpolationType,PressureInterpolationType, &
      & PressureInterpolationType],Err)
    CALL CMISSBasis_QuadratureNumberOfGaussXiSet(HydrostaticPressureBasis2,[NumberOfGaussXi,NumberOfGaussXi, &
      & NumberOfGaussXi],Err)
    !Finish the creation of the pressure basis
    CALL CMISSBasis_CreateFinish(HydrostaticPressureBasis2,Err)
  ENDIF
  
  !Start the creation of a generated mesh in the first region
  PRINT *, ' == >> CREATING GENERATED MESH(1) << == '
  CALL CMISSGeneratedMesh_Initialise(GeneratedMesh1,Err)
  !Set up a regular x*y*z mesh
  CALL CMISSGeneratedMesh_CreateStart(GeneratedMesh1UserNumber,Region1,GeneratedMesh1,Err)
  CALL CMISSGeneratedMesh_TypeSet(GeneratedMesh1,CMISS_GENERATED_MESH_REGULAR_MESH_TYPE,Err)
  IF(HydrostaticPressureBasis) THEN
    CALL CMISSGeneratedMesh_BasisSet(GeneratedMesh1,[DisplacementBasis1,HydrostaticPressureBasis1],Err)
  ELSE
    CALL CMISSGeneratedMesh_BasisSet(GeneratedMesh1,[DisplacementBasis1],Err)
  ENDIF
  !Define the mesh on the first region
  CALL CMISSGeneratedMesh_ExtentSet(GeneratedMesh1,[Width,Height,Length],Err) 
  CALL CMISSGeneratedMesh_OriginSet(GeneratedMesh1,[0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP],Err) 
  CALL CMISSGeneratedMesh_NumberOfElementsSet(GeneratedMesh1,[NumberGlobalXElements, &
    & NumberGlobalYElements,NumberGlobalZElements],Err)
  !Finish the creation of a generated mesh in the first region 
  CALL CMISSMesh_Initialise(Mesh1,Err) 
  CALL CMISSGeneratedMesh_CreateFinish(GeneratedMesh1,Mesh1UserNumber,Mesh1,Err)

  !Start the creation of a generated mesh in the second region
  PRINT *, ' == >> CREATING GENERATED MESH(2) << == '
  CALL CMISSGeneratedMesh_Initialise(GeneratedMesh2,Err)
  !Set up a regular x*y*z mesh
  CALL CMISSGeneratedMesh_CreateStart(GeneratedMesh2UserNumber,Region2,GeneratedMesh2,Err)
  CALL CMISSGeneratedMesh_TypeSet(GeneratedMesh2,CMISS_GENERATED_MESH_REGULAR_MESH_TYPE,Err)
  IF(HydrostaticPressureBasis) THEN
    CALL CMISSGeneratedMesh_BasisSet(GeneratedMesh2,[DisplacementBasis2,HydrostaticPressureBasis2],Err)
  ELSE
    CALL CMISSGeneratedMesh_BasisSet(GeneratedMesh2,[DisplacementBasis2],Err)
  ENDIF
  !Define the mesh on the first region
  IF(Uncoupled.EQV..TRUE.) THEN
    CALL CMISSGeneratedMesh_OriginSet(GeneratedMesh2,[0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP],Err) 
  ELSE
    CALL CMISSGeneratedMesh_OriginSet(GeneratedMesh2,[Width,0.0_CMISSDP,0.0_CMISSDP],Err) 
  ENDIF
  CALL CMISSGeneratedMesh_ExtentSet(GeneratedMesh2,[Width,Height,Length],Err) 
  CALL CMISSGeneratedMesh_NumberOfElementsSet(GeneratedMesh2,[NumberGlobalXElements, &
    & NumberGlobalYElements,NumberGlobalZElements],Err)
  !Finish the creation of a generated mesh in the first region 
  CALL CMISSMesh_Initialise(Mesh2,Err) 
  CALL CMISSGeneratedMesh_CreateFinish(GeneratedMesh2,Mesh2UserNumber,Mesh2,Err)

  !Create a decomposition for mesh1
  PRINT *, ' == >> CREATING MESH(1) DECOMPOSITION << == '
  CALL CMISSDecomposition_Initialise(Decomposition1,Err)
  CALL CMISSDecomposition_CreateStart(Decomposition1UserNumber,Mesh1,Decomposition1,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL CMISSDecomposition_TypeSet(Decomposition1,CMISS_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL CMISSDecomposition_NumberOfDomainsSet(Decomposition1,NumberOfComputationalNodes,Err)
  !Finish the decomposition
  CALL CMISSDecomposition_CreateFinish(Decomposition1,Err)

  !Create a decomposition for mesh2
  PRINT *, ' == >> CREATING MESH(2) DECOMPOSITION << == '
  CALL CMISSDecomposition_Initialise(Decomposition2,Err)
  CALL CMISSDecomposition_CreateStart(Decomposition2UserNumber,Mesh2,Decomposition2,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL CMISSDecomposition_TypeSet(Decomposition2,CMISS_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL CMISSDecomposition_NumberOfDomainsSet(Decomposition2,NumberOfComputationalNodes,Err)
  !Finish the decomposition
  CALL CMISSDecomposition_CreateFinish(Decomposition2,Err)
  
  !Start to create a default (geometric) field on the first region
  PRINT *, ' == >> CREATING MESH(1) GEOMETRIC FIELD << == '
  CALL CMISSField_Initialise(GeometricField1,Err)
  CALL CMISSField_CreateStart(GeometricField1UserNumber,Region1,GeometricField1,Err)
  !Set the decomposition to use
  CALL CMISSField_MeshDecompositionSet(GeometricField1,Decomposition1,Err)
  CALL CMISSField_TypeSet(GeometricField1,CMISS_FIELD_GEOMETRIC_TYPE,Err)
  CALL CMISSField_NumberOfVariablesSet(GeometricField1,GeometricField1NumberOfVariables,Err)
  CALL CMISSField_VariableLabelSet(GeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,"Geometry",Err)
  IF(DisplacementInterpolationType==CMISS_BASIS_CUBIC_HERMITE_INTERPOLATION)THEN
    CALL CMISSField_ScalingTypeSet(GeometricField1,CMISS_FIELD_ARITHMETIC_MEAN_SCALING,Err)
  ELSE
    CALL CMISSField_ScalingTypeSet(GeometricField1,CMISS_FIELD_UNIT_SCALING,Err)
  ENDIF
  CALL CMISSField_NumberOfComponentsSet(GeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,GeometricField1NumberOfComponents,Err)
  !Set the domain to be used by the field components.
  CALL CMISSField_ComponentMeshComponentSet(GeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,1,1,Err)
  CALL CMISSField_ComponentMeshComponentSet(GeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,2,1,Err)
  CALL CMISSField_ComponentMeshComponentSet(GeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,3,1,Err)
  !Finish creating the first field
  CALL CMISSField_CreateFinish(GeometricField1,Err)

  !Start to create a default (geometric) field on the second region
  PRINT *, ' == >> CREATING MESH(2) GEOMETRIC FIELD << == '
  CALL CMISSField_Initialise(GeometricField2,Err)
  CALL CMISSField_CreateStart(GeometricField2UserNumber,Region2,GeometricField2,Err)
  !Set the decomposition to use
  CALL CMISSField_MeshDecompositionSet(GeometricField2,Decomposition2,Err)
  CALL CMISSField_TypeSet(GeometricField2,CMISS_FIELD_GEOMETRIC_TYPE,Err)
  CALL CMISSField_NumberOfVariablesSet(GeometricField2,GeometricField2NumberOfVariables,Err)
  CALL CMISSField_VariableLabelSet(GeometricField2,CMISS_FIELD_U_VARIABLE_TYPE,"Geometry",Err)
  IF(DisplacementInterpolationType==CMISS_BASIS_CUBIC_HERMITE_INTERPOLATION)THEN
    CALL CMISSField_ScalingTypeSet(GeometricField2,CMISS_FIELD_ARITHMETIC_MEAN_SCALING,Err)
  ELSE
    CALL CMISSField_ScalingTypeSet(GeometricField2,CMISS_FIELD_UNIT_SCALING,Err)
  ENDIF
  CALL CMISSField_NumberOfComponentsSet(GeometricField2,CMISS_FIELD_U_VARIABLE_TYPE,GeometricField2NumberOfComponents,Err)
  !Set the domain to be used by the field components.
  CALL CMISSField_ComponentMeshComponentSet(GeometricField2,CMISS_FIELD_U_VARIABLE_TYPE,1,1,Err)
  CALL CMISSField_ComponentMeshComponentSet(GeometricField2,CMISS_FIELD_U_VARIABLE_TYPE,2,1,Err)
  CALL CMISSField_ComponentMeshComponentSet(GeometricField2,CMISS_FIELD_U_VARIABLE_TYPE,3,1,Err)
  !Finish creating the first field
  CALL CMISSField_CreateFinish(GeometricField2,Err)

  !Update the geometric field parameters for the first field
  PRINT *, ' == >> INITIALIZE GEOMETRIC FIELD(1) FROM GENERATED MESH(1) << == '
  CALL CMISSGeneratedMesh_GeometricParametersCalculate(GeneratedMesh1,GeometricField1,Err)

  !Update the geometric field parameters for the second field
  PRINT *, ' == >> INITIALIZE GEOMETRIC FIELD(2) FROM GENERATED MESH(2) << == '
  CALL CMISSGeneratedMesh_GeometricParametersCalculate(GeneratedMesh2,GeometricField2,Err)

  !######################################### FINITE ELASTICITY ##########################################

  !Start to create a default (fibre) field on the first region
  PRINT *, ' == >> CREATING MESH(1) FIBRE FIELD << == '
  CALL CMISSField_Initialise(FibreField1,Err)
  CALL CMISSField_CreateStart(FibreField1UserNumber,Region1,FibreField1,Err)
  CALL CMISSField_TypeSet(FibreField1,CMISS_FIELD_FIBRE_TYPE,Err)
  !Set the decomposition to use
  CALL CMISSField_MeshDecompositionSet(FibreField1,Decomposition1,Err)
  CALL CMISSField_GeometricFieldSet(FibreField1,GeometricField1,Err)
  CALL CMISSField_NumberOfVariablesSet(FibreField1,FibreField1NumberOfVariables,Err)
  CALL CMISSField_VariableLabelSet(FibreField1,CMISS_FIELD_U_VARIABLE_TYPE,"Fibre",Err)
  IF(DisplacementInterpolationType==CMISS_BASIS_CUBIC_HERMITE_INTERPOLATION)THEN
    CALL CMISSField_ScalingTypeSet(FibreField1,CMISS_FIELD_ARITHMETIC_MEAN_SCALING,Err)
  ELSE
    CALL CMISSField_ScalingTypeSet(FibreField1,CMISS_FIELD_UNIT_SCALING,Err)
  ENDIF
  CALL CMISSField_NumberOfComponentsSet(FibreField1,CMISS_FIELD_U_VARIABLE_TYPE,FibreField1NumberOfComponents,Err)  
  !Set the domain to be used by the field components.
  CALL CMISSField_ComponentMeshComponentSet(FibreField1,CMISS_FIELD_U_VARIABLE_TYPE,1,1,Err)
  CALL CMISSField_ComponentMeshComponentSet(FibreField1,CMISS_FIELD_U_VARIABLE_TYPE,2,1,Err)
  CALL CMISSField_ComponentMeshComponentSet(FibreField1,CMISS_FIELD_U_VARIABLE_TYPE,3,1,Err)
  !Finish creating the first field
  CALL CMISSField_CreateFinish(FibreField1,Err)

  !Start to create a default (fibre) field on the second region
  PRINT *, ' == >> CREATING MESH(2) FIBRE FIELD << == '
  CALL CMISSField_Initialise(FibreField2,Err)
  CALL CMISSField_CreateStart(FibreField2UserNumber,Region2,FibreField2,Err)
  CALL CMISSField_TypeSet(FibreField2,CMISS_FIELD_FIBRE_TYPE,Err)
  !Set the decomposition to use
  CALL CMISSField_MeshDecompositionSet(FibreField2,Decomposition2,Err)
  CALL CMISSField_GeometricFieldSet(FibreField2,GeometricField2,Err)
  CALL CMISSField_NumberOfVariablesSet(FibreField2,FibreField2NumberOfVariables,Err)
  CALL CMISSField_VariableLabelSet(FibreField2,CMISS_FIELD_U_VARIABLE_TYPE,"Fibre",Err)
  IF(DisplacementInterpolationType==CMISS_BASIS_CUBIC_HERMITE_INTERPOLATION)THEN
    CALL CMISSField_ScalingTypeSet(FibreField2,CMISS_FIELD_ARITHMETIC_MEAN_SCALING,Err)
  ELSE
    CALL CMISSField_ScalingTypeSet(FibreField2,CMISS_FIELD_UNIT_SCALING,Err)
  ENDIF
  CALL CMISSField_NumberOfComponentsSet(FibreField2,CMISS_FIELD_U_VARIABLE_TYPE,FibreField2NumberOfComponents,Err)  
  !Set the domain to be used by the field components.
  CALL CMISSField_ComponentMeshComponentSet(FibreField2,CMISS_FIELD_U_VARIABLE_TYPE,1,1,Err)
  CALL CMISSField_ComponentMeshComponentSet(FibreField2,CMISS_FIELD_U_VARIABLE_TYPE,2,1,Err)
  CALL CMISSField_ComponentMeshComponentSet(FibreField2,CMISS_FIELD_U_VARIABLE_TYPE,3,1,Err)
  !Finish creating the first field
  CALL CMISSField_CreateFinish(FibreField2,Err)

  !Start to create a material field on the first region
  PRINT *, ' == >> CREATING MESH(1) MATERIAL FIELD << == '
  CALL CMISSField_Initialise(MaterialField1,Err)
  CALL CMISSField_CreateStart(MaterialField1UserNumber,Region1,MaterialField1,Err)
  CALL CMISSField_TypeSet(MaterialField1,CMISS_FIELD_MATERIAL_TYPE,Err)
  !Set the decomposition to use
  CALL CMISSField_MeshDecompositionSet(MaterialField1,Decomposition1,Err)
  CALL CMISSField_GeometricFieldSet(MaterialField1,GeometricField1,Err)
  CALL CMISSField_NumberOfVariablesSet(MaterialField1,MaterialField1NumberOfVariables,Err)
  CALL CMISSField_VariableLabelSet(MaterialField1,CMISS_FIELD_U_VARIABLE_TYPE,"Material1",Err)
  IF(DisplacementInterpolationType==CMISS_BASIS_CUBIC_HERMITE_INTERPOLATION)THEN
    CALL CMISSField_ScalingTypeSet(MaterialField1,CMISS_FIELD_ARITHMETIC_MEAN_SCALING,Err)
  ELSE
    CALL CMISSField_ScalingTypeSet(MaterialField1,CMISS_FIELD_UNIT_SCALING,Err)
  ENDIF
  CALL CMISSField_NumberOfComponentsSet(MaterialField1,CMISS_FIELD_U_VARIABLE_TYPE,MaterialField1NumberOfComponents,Err)  
  !Set the domain to be used by the field components.
  CALL CMISSField_ComponentMeshComponentSet(MaterialField1,CMISS_FIELD_U_VARIABLE_TYPE,1,1,Err)
  CALL CMISSField_ComponentMeshComponentSet(MaterialField1,CMISS_FIELD_U_VARIABLE_TYPE,2,1,Err)
  CALL CMISSField_ComponentInterpolationSet(MaterialField1,CMISS_FIELD_U_VARIABLE_TYPE,1, &
    & CMISS_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  CALL CMISSField_ComponentInterpolationSet(MaterialField1,CMISS_FIELD_U_VARIABLE_TYPE,2, &
    & CMISS_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  IF (Incompressible1.eqv..FALSE.) THEN
    CALL CMISSField_ComponentMeshComponentSet(MaterialField1,CMISS_FIELD_U_VARIABLE_TYPE,3,1,Err)
    CALL CMISSField_ComponentInterpolationSet(MaterialField1,CMISS_FIELD_U_VARIABLE_TYPE,3, &
      & CMISS_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  ENDIF
  !Finish creating the first field
  CALL CMISSField_CreateFinish(MaterialField1,Err)

  !Start to create a material field on the second region
  PRINT *, ' == >> CREATING MESH(2) MATERIAL FIELD << == '
  CALL CMISSField_Initialise(MaterialField2,Err)
  CALL CMISSField_CreateStart(MaterialField2UserNumber,Region2,MaterialField2,Err)
  CALL CMISSField_TypeSet(MaterialField2,CMISS_FIELD_MATERIAL_TYPE,Err)
  !Set the decomposition to use
  CALL CMISSField_MeshDecompositionSet(MaterialField2,Decomposition2,Err)
  CALL CMISSField_GeometricFieldSet(MaterialField2,GeometricField2,Err)
  CALL CMISSField_NumberOfVariablesSet(MaterialField2,MaterialField2NumberOfVariables,Err)
  CALL CMISSField_VariableLabelSet(MaterialField2,CMISS_FIELD_U_VARIABLE_TYPE,"Material2",Err)
  IF(DisplacementInterpolationType==CMISS_BASIS_CUBIC_HERMITE_INTERPOLATION)THEN
    CALL CMISSField_ScalingTypeSet(MaterialField2,CMISS_FIELD_ARITHMETIC_MEAN_SCALING,Err)
  ELSE
    CALL CMISSField_ScalingTypeSet(MaterialField2,CMISS_FIELD_UNIT_SCALING,Err)
  ENDIF
  CALL CMISSField_NumberOfComponentsSet(MaterialField2,CMISS_FIELD_U_VARIABLE_TYPE,MaterialField2NumberOfComponents,Err)  
  !Set the domain to be used by the field components.
  CALL CMISSField_ComponentMeshComponentSet(MaterialField2,CMISS_FIELD_U_VARIABLE_TYPE,1,1,Err)
  CALL CMISSField_ComponentMeshComponentSet(MaterialField2,CMISS_FIELD_U_VARIABLE_TYPE,2,1,Err)
  CALL CMISSField_ComponentInterpolationSet(MaterialField2,CMISS_FIELD_U_VARIABLE_TYPE,1, &
    & CMISS_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  CALL CMISSField_ComponentInterpolationSet(MaterialField2,CMISS_FIELD_U_VARIABLE_TYPE,2, &
    & CMISS_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  IF (Incompressible2.eqv..FALSE.) THEN
    CALL CMISSField_ComponentMeshComponentSet(MaterialField2,CMISS_FIELD_U_VARIABLE_TYPE,3,1,Err)
    CALL CMISSField_ComponentInterpolationSet(MaterialField2,CMISS_FIELD_U_VARIABLE_TYPE,3, &
      & CMISS_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  ENDIF
  !Finish creating the first field
  CALL CMISSField_CreateFinish(MaterialField2,Err)

  !Set Mooney-Rivlin constants c10 and c01
  PRINT *, ' == >> INITIALIZE MATERIAL FIELD(1) << == '
  CALL CMISSField_ComponentValuesInitialise(MaterialField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,2.0_CMISSDP,Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 2,6.0_CMISSDP,Err)
  IF (Incompressible1.eqv..FALSE.) THEN
    CALL CMISSField_ComponentValuesInitialise(MaterialField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
      & 3,1.0e3_CMISSDP,Err)
  ENDIF

  !Set Mooney-Rivlin constants c10 and c01
  PRINT *, ' == >> INITIALIZE MATERIAL FIELD(2) << == '
  CALL CMISSField_ComponentValuesInitialise(MaterialField2,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,2.0_CMISSDP,Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialField2,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 2,6.0_CMISSDP,Err)
  IF (Incompressible2.eqv..FALSE.) THEN
    CALL CMISSField_ComponentValuesInitialise(MaterialField2,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
      & 3,1.0e3_CMISSDP,Err)
  ENDIF

  !Start to create a dependent field on the first region
  PRINT *, ' == >> CREATING MESH(1) DEPENDENT FIELD << == '
  CALL CMISSField_Initialise(DependentField1,Err)
  CALL CMISSField_CreateStart(DependentField1UserNumber,Region1,DependentField1,Err)
  CALL CMISSField_TypeSet(DependentField1,CMISS_FIELD_GENERAL_TYPE,Err)
  !Set the decomposition to use
  CALL CMISSField_MeshDecompositionSet(DependentField1,Decomposition1,Err)
  CALL CMISSField_GeometricFieldSet(DependentField1,GeometricField1,Err)
  CALL CMISSField_DependentTypeSet(DependentField1,CMISS_FIELD_DEPENDENT_TYPE,Err) 
  CALL CMISSField_NumberOfVariablesSet(DependentField1,DependentField1NumberOfVariables,Err)
  CALL CMISSField_VariableLabelSet(DependentField1,CMISS_FIELD_U_VARIABLE_TYPE,"Dependent1",Err)
  CALL CMISSField_VariableLabelSet(DependentField1,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,"Dependent1DelUDelN",Err)
  IF(DisplacementInterpolationType==CMISS_BASIS_CUBIC_HERMITE_INTERPOLATION)THEN
    CALL CMISSField_ScalingTypeSet(DependentField1,CMISS_FIELD_ARITHMETIC_MEAN_SCALING,Err)
  ELSE
    CALL CMISSField_ScalingTypeSet(DependentField1,CMISS_FIELD_UNIT_SCALING,Err)
  ENDIF
  CALL CMISSField_NumberOfComponentsSet(DependentField1,CMISS_FIELD_U_VARIABLE_TYPE,DependentField1NumberOfComponents,Err)  
  CALL CMISSField_NumberOfComponentsSet(DependentField1,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,DependentField1NumberOfComponents,Err)
  !Set the domain to be used by the field components.
  CALL CMISSField_ComponentMeshComponentSet(DependentField1,CMISS_FIELD_U_VARIABLE_TYPE,1,1,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField1,CMISS_FIELD_U_VARIABLE_TYPE,2,1,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField1,CMISS_FIELD_U_VARIABLE_TYPE,3,1,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField1,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,1,1,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField1,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,2,1,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField1,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,3,1,Err)
  IF (Incompressible1.eqv..TRUE.) THEN
    IF(HydrostaticPressureBasis) THEN
      CALL CMISSField_ComponentMeshComponentSet(DependentField1,CMISS_FIELD_U_VARIABLE_TYPE,4,2,Err)
      CALL CMISSField_ComponentMeshComponentSet(DependentField1,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,4,2,Err)
      !Set the pressure to be nodally based and use the second mesh component if required
      CALL CMISSField_ComponentInterpolationSet(DependentField1,CMISS_FIELD_U_VARIABLE_TYPE,4, &
        & CMISS_FIELD_NODE_BASED_INTERPOLATION,Err)
      CALL CMISSField_ComponentInterpolationSet(DependentField1,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,4, &
        & CMISS_FIELD_NODE_BASED_INTERPOLATION,Err)
    ELSE
      CALL CMISSField_ComponentMeshComponentSet(DependentField1,CMISS_FIELD_U_VARIABLE_TYPE,4,1,Err)
      CALL CMISSField_ComponentMeshComponentSet(DependentField1,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,4,1,Err)
      CALL CMISSField_ComponentInterpolationSet(DependentField1,CMISS_FIELD_U_VARIABLE_TYPE,4, &
        & CMISS_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
      CALL CMISSField_ComponentInterpolationSet(DependentField1,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,4, &
        & CMISS_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
    ENDIF
  ENDIF
  !Finish creating the first field
  CALL CMISSField_CreateFinish(DependentField1,Err)

  !Start to create a dependent field on the second region
  PRINT *, ' == >> CREATING MESH(2) DEPENDENT FIELD << == '
  CALL CMISSField_Initialise(DependentField2,Err)
  CALL CMISSField_CreateStart(DependentField2UserNumber,Region2,DependentField2,Err)
  CALL CMISSField_TypeSet(DependentField2,CMISS_FIELD_GENERAL_TYPE,Err)
  !Set the decomposition to use
  CALL CMISSField_MeshDecompositionSet(DependentField2,Decomposition2,Err)
  CALL CMISSField_GeometricFieldSet(DependentField2,GeometricField2,Err)
  CALL CMISSField_DependentTypeSet(DependentField2,CMISS_FIELD_DEPENDENT_TYPE,Err) 
  CALL CMISSField_NumberOfVariablesSet(DependentField2,DependentField2NumberOfVariables,Err)
  CALL CMISSField_VariableLabelSet(DependentField2,CMISS_FIELD_U_VARIABLE_TYPE,"Dependent2",Err)
  CALL CMISSField_VariableLabelSet(DependentField2,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,"Dependent2DelUDelN",Err)
  IF(DisplacementInterpolationType==CMISS_BASIS_CUBIC_HERMITE_INTERPOLATION)THEN
    CALL CMISSField_ScalingTypeSet(DependentField2,CMISS_FIELD_ARITHMETIC_MEAN_SCALING,Err)
  ELSE
    CALL CMISSField_ScalingTypeSet(DependentField2,CMISS_FIELD_UNIT_SCALING,Err)
  ENDIF
  CALL CMISSField_NumberOfComponentsSet(DependentField2,CMISS_FIELD_U_VARIABLE_TYPE,DependentField2NumberOfComponents,Err)  
  CALL CMISSField_NumberOfComponentsSet(DependentField2,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,DependentField2NumberOfComponents,Err)
  !Set the domain to be used by the field components.
  CALL CMISSField_ComponentMeshComponentSet(DependentField2,CMISS_FIELD_U_VARIABLE_TYPE,1,1,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField2,CMISS_FIELD_U_VARIABLE_TYPE,2,1,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField2,CMISS_FIELD_U_VARIABLE_TYPE,3,1,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField2,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,1,1,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField2,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,2,1,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField2,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,3,1,Err)
  IF (Incompressible2.eqv..TRUE.) THEN
    IF(HydrostaticPressureBasis) THEN
      CALL CMISSField_ComponentMeshComponentSet(DependentField2,CMISS_FIELD_U_VARIABLE_TYPE,4,2,Err)
      CALL CMISSField_ComponentMeshComponentSet(DependentField2,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,4,2,Err)
      !Set the pressure to be nodally based and use the second mesh component if required
      CALL CMISSField_ComponentInterpolationSet(DependentField2,CMISS_FIELD_U_VARIABLE_TYPE,4, &
        & CMISS_FIELD_NODE_BASED_INTERPOLATION,Err)
      CALL CMISSField_ComponentInterpolationSet(DependentField2,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,4, &
        & CMISS_FIELD_NODE_BASED_INTERPOLATION,Err)
    ELSE
      CALL CMISSField_ComponentMeshComponentSet(DependentField2,CMISS_FIELD_U_VARIABLE_TYPE,4,1,Err)
      CALL CMISSField_ComponentMeshComponentSet(DependentField2,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,4,1,Err)
      CALL CMISSField_ComponentInterpolationSet(DependentField2,CMISS_FIELD_U_VARIABLE_TYPE,4, &
        & CMISS_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
      CALL CMISSField_ComponentInterpolationSet(DependentField2,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,4, &
        & CMISS_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
    ENDIF
  ENDIF
  !Finish creating the first field
  CALL CMISSField_CreateFinish(DependentField2,Err)

  !Initialise dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
  PRINT *, ' == >> INITIALIZE DEPENDENT FIELD(1) << == '
  CALL CMISSField_ParametersToFieldParametersComponentCopy(GeometricField1,CMISS_FIELD_U_VARIABLE_TYPE, &
    & CMISS_FIELD_VALUES_SET_TYPE,1,DependentField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,Err)
  CALL CMISSField_ParametersToFieldParametersComponentCopy(GeometricField1,CMISS_FIELD_U_VARIABLE_TYPE, &
    & CMISS_FIELD_VALUES_SET_TYPE,2,DependentField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2,Err)
  CALL CMISSField_ParametersToFieldParametersComponentCopy(GeometricField1,CMISS_FIELD_U_VARIABLE_TYPE, &
    & CMISS_FIELD_VALUES_SET_TYPE,3,DependentField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,3,Err)
  IF (Incompressible1.eqv..TRUE.) THEN
    CALL CMISSField_ComponentValuesInitialise(DependentField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
      & 4,-8.0_CMISSDP,Err)
  ENDIF

  !Initialise dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
  PRINT *, ' == >> INITIALIZE DEPENDENT FIELD(2) << == '
  CALL CMISSField_ParametersToFieldParametersComponentCopy(GeometricField2,CMISS_FIELD_U_VARIABLE_TYPE, &
    & CMISS_FIELD_VALUES_SET_TYPE,1,DependentField2,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,Err)
  CALL CMISSField_ParametersToFieldParametersComponentCopy(GeometricField2,CMISS_FIELD_U_VARIABLE_TYPE, &
    & CMISS_FIELD_VALUES_SET_TYPE,2,DependentField2,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2,Err)
  CALL CMISSField_ParametersToFieldParametersComponentCopy(GeometricField2,CMISS_FIELD_U_VARIABLE_TYPE, &
    & CMISS_FIELD_VALUES_SET_TYPE,3,DependentField2,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,3,Err)
  IF (Incompressible2.eqv..TRUE.) THEN
    CALL CMISSField_ComponentValuesInitialise(DependentField2,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
      & 4,-8.0_CMISSDP,Err)
  ENDIF

  !Create the equations set for the first region
  PRINT *, ' == >> CREATING EQUATION SET(1) << == '
  CALL CMISSField_Initialise(EquationsSetField1,Err)
  CALL CMISSEquationsSet_Initialise(EquationsSet1,Err)
  !Set the equations set to be a Finite Elasticity problem
  IF (Incompressible1.eqv..TRUE.) THEN
    CALL CMISSEquationsSet_CreateStart(EquationsSet1UserNumber,Region1,FibreField1,CMISS_EQUATIONS_SET_ELASTICITY_CLASS, &
      & CMISS_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,CMISS_EQUATIONS_SET_NO_SUBTYPE,EquationsSetField1UserNumber,&
      & EquationsSetField1,EquationsSet1,Err)
  ELSE
    CALL CMISSEquationsSet_CreateStart(EquationsSet1UserNumber,Region1,FibreField1,CMISS_EQUATIONS_SET_ELASTICITY_CLASS, &
      & CMISS_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,CMISS_EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE, &
      & EquationsSetField1UserNumber,EquationsSetField1,EquationsSet1,Err)
  ENDIF

  !Finish creating the equations set
  CALL CMISSEquationsSet_CreateFinish(EquationsSet1,Err)

  !Create the equations set for the second region
  PRINT *, ' == >> CREATING EQUATION SET(2) << == '
  CALL CMISSField_Initialise(EquationsSetField2,Err)
  CALL CMISSEquationsSet_Initialise(EquationsSet2,Err)
  !Set the equations set to be a Finite Elasticity problem
  IF (Incompressible2.eqv..TRUE.) THEN
    CALL CMISSEquationsSet_CreateStart(EquationsSet2UserNumber,Region2,FibreField2,CMISS_EQUATIONS_SET_ELASTICITY_CLASS, &
      & CMISS_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,CMISS_EQUATIONS_SET_NO_SUBTYPE,EquationsSetField2UserNumber,&
      & EquationsSetField2,EquationsSet2,Err)
  ELSE
    CALL CMISSEquationsSet_CreateStart(EquationsSet2UserNumber,Region2,FibreField2,CMISS_EQUATIONS_SET_ELASTICITY_CLASS, &
      & CMISS_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,CMISS_EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE, &
      & EquationsSetField2UserNumber,EquationsSetField2,EquationsSet2,Err)
  ENDIF
  !Finish creating the equations set
  CALL CMISSEquationsSet_CreateFinish(EquationsSet2,Err)

  !Create the equations set dependent field variables for the first equations set
  PRINT *, ' == >> CREATING DEPENDENT FIELD EQUATION SET(1) << == '
  CALL CMISSEquationsSet_DependentCreateStart(EquationsSet1,DependentField1UserNumber,DependentField1,Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSet_DependentCreateFinish(EquationsSet1,Err)

  !Create the equations set dependent field variables for the second equations set
  PRINT *, ' == >> CREATING DEPENDENT FIELD EQUATION SET(2) << == '
  CALL CMISSEquationsSet_DependentCreateStart(EquationsSet2,DependentField2UserNumber,DependentField2,Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSet_DependentCreateFinish(EquationsSet2,Err)

  !Create the equations set dependent field variables for the first equations set
  PRINT *, ' == >> CREATING MATERIAL FIELD EQUATION SET(1) << == '
  CALL CMISSEquationsSet_MaterialsCreateStart(EquationsSet1,MaterialField1UserNumber,MaterialField1,Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSet_MaterialsCreateFinish(EquationsSet1,Err)

  !Create the equations set dependent field variables for the second equations set
  PRINT *, ' == >> CREATING MATERIAL FIELD EQUATION SET(2) << == '
  CALL CMISSEquationsSet_MaterialsCreateStart(EquationsSet2,MaterialField2UserNumber,MaterialField2,Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSet_MaterialsCreateFinish(EquationsSet2,Err)

  !Create the equations set equations for the first equations set
  PRINT *, ' == >> CREATING EQUATIONS(1) << == '
  CALL CMISSEquations_Initialise(Equations1,Err)
  CALL CMISSEquationsSet_EquationsCreateStart(EquationsSet1,Equations1,Err)
  !Set the equations matrices sparsity type
  CALL CMISSEquations_SparsityTypeSet(Equations1,CMISS_EQUATIONS_SPARSE_MATRICES,Err)
  !CALL CMISSEquations_SparsityTypeSet(Equations1,CMISS_EQUATIONS_FULL_MATRICES,Err)
  !Set the equations set output
  CALL CMISSEquations_OutputTypeSet(Equations1,CMISS_EQUATIONS_NO_OUTPUT,Err)
  !CALL CMISSEquations_OutputTypeSet(Equations1,CMISS_EQUATIONS_TIMING_OUTPUT,Err)
  !CALL CMISSEquations_OutputTypeSet(Equations1,CMISS_EQUATIONS_MATRIX_OUTPUT,Err)
  !CALL CMISSEquations_OutputTypeSet(Equations1,CMISS_EQUATIONS_ELEMENT_MATRIX_OUTPUT,Err)
  !Finish the equations set equations
  CALL CMISSEquationsSet_EquationsCreateFinish(EquationsSet1,Err)

  !Create the equations set equations for the second equations set
  PRINT *, ' == >> CREATING EQUATIONS(2) << == '
  CALL CMISSEquations_Initialise(Equations2,Err)
  CALL CMISSEquationsSet_EquationsCreateStart(EquationsSet2,Equations2,Err)
  !Set the equations matrices sparsity type
  CALL CMISSEquations_SparsityTypeSet(Equations2,CMISS_EQUATIONS_SPARSE_MATRICES,Err)
  !CALL CMISSEquations_SparsityTypeSet(Equations2,CMISS_EQUATIONS_FULL_MATRICES,Err)
  !Set the equations set output
  CALL CMISSEquations_OutputTypeSet(Equations2,CMISS_EQUATIONS_NO_OUTPUT,Err)
  !CALL CMISSEquations_OutputTypeSet(Equations2,CMISS_EQUATIONS_TIMING_OUTPUT,Err)
  !CALL CMISSEquations_OutputTypeSet(Equations2,CMISS_EQUATIONS_MATRIX_OUTPUT,Err)
  !CALL CMISSEquations_OutputTypeSet(Equations2,CMISS_EQUATIONS_ELEMENT_MATRIX_OUTPUT,Err)
  !Finish the equations set equations
  CALL CMISSEquationsSet_EquationsCreateFinish(EquationsSet2,Err)

  !######################################### START CREATING INTERFACE ##########################################

  !Create an interface between the two meshes
  PRINT *, ' == >> CREATING INTERFACE << == '
  CALL CMISSInterface_Initialise(Interface,Err)
  CALL CMISSInterface_CreateStart(InterfaceUserNumber,WorldRegion,Interface,Err)
  CALL CMISSInterface_LabelSet(Interface,"Interface",Err)
  !Add in the two meshes
  CALL CMISSInterface_MeshAdd(Interface,Mesh1,Mesh1Index,Err)
  CALL CMISSInterface_MeshAdd(Interface,Mesh2,Mesh2Index,Err)
  !Finish creating the interface
  CALL CMISSInterface_CreateFinish(Interface,Err)

  PRINT *, ' == >> CREATING INTERFACE BASIS << == '
  CALL CMISSBasis_Initialise(InterfaceBasis,Err)
  CALL CMISSBasis_CreateStart(InterfaceBasisUserNumber,InterfaceBasis,Err)
  CALL CMISSBasis_NumberOfXiSet(InterfaceBasis,2,Err)
  CALL CMISSBasis_InterpolationXiSet(InterfaceBasis,[DisplacementInterpolationType, &
    & DisplacementInterpolationType],Err)
  !Finish the creation of the basis
  CALL CMISSBasis_CreateFinish(InterfaceBasis,Err)

  PRINT *, ' == >> CREATING INTERFACE MAPPING BASIS << == '
  CALL CMISSBasis_Initialise(InterfaceMappingBasis,Err)
  CALL CMISSBasis_CreateStart(InterfaceMappingBasisUserNumber,InterfaceMappingBasis,Err)
  CALL CMISSBasis_NumberOfXiSet(InterfaceMappingBasis,2,Err)
  CALL CMISSBasis_InterpolationXiSet(InterfaceMappingBasis,[CMISS_BASIS_LINEAR_LAGRANGE_INTERPOLATION, &
    & CMISS_BASIS_LINEAR_LAGRANGE_INTERPOLATION],Err)
  !Finish the creation of the basis
  CALL CMISSBasis_CreateFinish(InterfaceMappingBasis,Err)
  
  !Start the creation of a generated mesh for the interface
  PRINT *, ' == >> CREATING INTERFACE GENERATED MESH << == '
  CALL CMISSGeneratedMesh_Initialise(InterfaceGeneratedMesh,Err)
  CALL CMISSGeneratedMesh_CreateStart(InterfaceGeneratedMeshUserNumber,Interface,InterfaceGeneratedMesh,Err)
  !Set up a regular x*y*z mesh
  CALL CMISSGeneratedMesh_TypeSet(InterfaceGeneratedMesh,CMISS_GENERATED_MESH_REGULAR_MESH_TYPE,Err)
  !Set the default basis
  CALL CMISSGeneratedMesh_BasisSet(InterfaceGeneratedMesh,InterfaceBasis,Err)
  CALL CMISSGeneratedMesh_OriginSet(InterfaceGeneratedMesh,[Width,0.0_CMISSDP,0.0_CMISSDP],Err)
  CALL CMISSGeneratedMesh_ExtentSet(InterfaceGeneratedMesh,[0.0_CMISSDP,Height,Length],Err)
  CALL CMISSGeneratedMesh_NumberOfElementsSet(InterfaceGeneratedMesh,[NumberGlobalYElements, &
    & NumberGlobalZElements],Err)
  !Finish the creation of a generated mesh in interface
  CALL CMISSMesh_Initialise(InterfaceMesh,Err)
  CALL CMISSGeneratedMesh_CreateFinish(InterfaceGeneratedMesh,InterfaceMeshUserNumber,InterfaceMesh,Err)

  !Couple the interface meshes
  PRINT *, ' == >> CREATING INTERFACE MESHES CONNECTIVITY << == '
  CALL CMISSInterfaceMeshConnectivity_Initialise(InterfaceMeshConnectivity,Err)
  CALL CMISSInterfaceMeshConnectivity_CreateStart(Interface,InterfaceMesh,InterfaceMeshConnectivity,Err)
  CALL CMISSInterfaceMeshConnectivity_BasisSet(InterfaceMeshConnectivity,InterfaceMappingBasis,Err)
  SELECT CASE(DisplacementInterpolationType)
  CASE(1,4)
    NumberOfNodeXi=2
  CASE(2)
    NumberOfNodeXi=3
  CASE(3)
    NumberOfNodeXi=4
  CASE DEFAULT
    CALL HANDLE_ERROR("Invalid interpolation type.")
  END SELECT
  DO y_element_idx=1,NumberGlobalYElements
    DO z_element_idx=1,NumberGlobalZElements
      !Map the interface element to the elements in mesh 1
      CALL CMISSInterfaceMeshConnectivity_ElementNumberSet(InterfaceMeshConnectivity, &
        & y_element_idx+(z_element_idx-1)*NumberGlobalYElements,Mesh1Index, &
        y_element_idx*NumberGlobalXElements+(z_element_idx-1)*NumberGlobalXElements* &
        & NumberGlobalYElements,Err)
      XI3 = [ 1.0_CMISSDP, 0.0_CMISSDP, 0.0_CMISSDP ]
      CALL CMISSInterfaceMeshConnectivity_ElementXiSet(InterfaceMeshConnectivity,y_element_idx+ &
        & (z_element_idx-1)*NumberGlobalYElements,Mesh1Index,y_element_idx* &
        & NumberGlobalXElements+(z_element_idx-1)*NumberGlobalXElements* &
        & NumberGlobalYElements,1,1,XI3,Err)
      XI3 = [ 1.0_CMISSDP, 1.0_CMISSDP, 0.0_CMISSDP ]
      CALL CMISSInterfaceMeshConnectivity_ElementXiSet(InterfaceMeshConnectivity,y_element_idx+ &
        & (z_element_idx-1)*NumberGlobalYElements,Mesh1Index,y_element_idx* &
        & NumberGlobalXElements+(z_element_idx-1)*NumberGlobalXElements* &
        & NumberGlobalYElements,2,1,XI3,Err)
      XI3 = [ 1.0_CMISSDP, 0.0_CMISSDP, 1.0_CMISSDP ]
      CALL CMISSInterfaceMeshConnectivity_ElementXiSet(InterfaceMeshConnectivity,y_element_idx+ &
        & (z_element_idx-1)*NumberGlobalYElements,Mesh1Index,y_element_idx* &
        & NumberGlobalXElements+(z_element_idx-1)*NumberGlobalXElements* &
        & NumberGlobalYElements,3,1,XI3,Err)
      XI3 = [ 1.0_CMISSDP, 1.0_CMISSDP, 1.0_CMISSDP ]
      CALL CMISSInterfaceMeshConnectivity_ElementXiSet(InterfaceMeshConnectivity,y_element_idx+ &
        & (z_element_idx-1)*NumberGlobalYElements,Mesh1Index,y_element_idx* &
        & NumberGlobalXElements+(z_element_idx-1)*NumberGlobalXElements* &
        & NumberGlobalYElements,4,1,XI3,Err)
      !Map the interface element to the elements in mesh 2
      CALL CMISSInterfaceMeshConnectivity_ElementNumberSet(InterfaceMeshConnectivity,y_element_idx+ &
        & (z_element_idx-1)*NumberGlobalYElements,Mesh2Index,1+(y_element_idx-1)* &
        & NumberGlobalXElements+(z_element_idx-1)*NumberGlobalXElements* &
        & NumberGlobalYElements,Err)
      DO mesh_local_y_node = 1,NumberOfNodeXi-1
        DO mesh_local_z_node = 1,NumberOfNodeXi-1
          XI3 = [ 0.0_CMISSDP, &
            & REAL(mesh_local_y_node-1,CMISSDP)/REAL(NumberOfNodeXi-1,CMISSDP), &
            & REAL(mesh_local_z_node-1,CMISSDP)/REAL(NumberOfNodeXi-1,CMISSDP) ]
          CALL CMISSInterfaceMeshConnectivity_ElementXiSet(InterfaceMeshConnectivity,y_element_idx+ &
            & (z_element_idx-1)*NumberGlobalYElements,Mesh2Index,1+(y_element_idx-1)* &
            & NumberGlobalXElements+(z_element_idx-1)*NumberGlobalXElements* &
            & NumberGlobalYElements,1,1,XI3,Err)
          XI3 = [ 0.0_CMISSDP, &
            & REAL(mesh_local_y_node,CMISSDP)/REAL(NumberOfNodeXi-1,CMISSDP), &
            & REAL(mesh_local_z_node-1,CMISSDP)/REAL(NumberOfNodeXi-1,CMISSDP) ]
          CALL CMISSInterfaceMeshConnectivity_ElementXiSet(InterfaceMeshConnectivity,y_element_idx+ &
            & (z_element_idx-1)*NumberGlobalYElements,Mesh2Index,1+(y_element_idx-1)* &
            & NumberGlobalXElements+(z_element_idx-1)*NumberGlobalXElements* &
            & NumberGlobalYElements,2,1,XI3,Err)
          XI3 = [ 0.0_CMISSDP, &
            & REAL(mesh_local_y_node-1,CMISSDP)/REAL(NumberOfNodeXi-1,CMISSDP), &
            & REAL(mesh_local_z_node,CMISSDP)/REAL(NumberOfNodeXi-1,CMISSDP) ]
          CALL CMISSInterfaceMeshConnectivity_ElementXiSet(InterfaceMeshConnectivity,y_element_idx+ &
            & (z_element_idx-1)*NumberGlobalYElements,Mesh2Index,1+(y_element_idx-1)* &
            & NumberGlobalXElements+(z_element_idx-1)*NumberGlobalXElements* &
            & NumberGlobalYElements,3,1,XI3,Err)
          XI3 = [ 0.0_CMISSDP, &
            & REAL(mesh_local_y_node,CMISSDP)/REAL(NumberOfNodeXi-1,CMISSDP), &
            & REAL(mesh_local_z_node,CMISSDP)/REAL(NumberOfNodeXi-1,CMISSDP) ]
          CALL CMISSInterfaceMeshConnectivity_ElementXiSet(InterfaceMeshConnectivity,y_element_idx+ &
            & (z_element_idx-1)*NumberGlobalYElements,Mesh2Index,1+(y_element_idx-1)* &
            & NumberGlobalXElements+(z_element_idx-1)*NumberGlobalXElements* &
            & NumberGlobalYElements,4,1,XI3,Err)
        ENDDO !mesh_local_z_node
      ENDDO !mesh_local_y_node
    ENDDO !z_element_idx
  ENDDO !y_element_idx
  CALL CMISSInterfaceMeshConnectivity_CreateFinish(InterfaceMeshConnectivity,Err)

  !Create a decomposition for the interface mesh
  PRINT *, ' == >> CREATING INTERFACE DECOMPOSITION << == '
  CALL CMISSDecomposition_Initialise(InterfaceDecomposition,Err)
  CALL CMISSDecomposition_CreateStart(InterfaceDecompositionUserNumber,InterfaceMesh,InterfaceDecomposition,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL CMISSDecomposition_TypeSet(InterfaceDecomposition,CMISS_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL CMISSDecomposition_NumberOfDomainsSet(InterfaceDecomposition,NumberOfComputationalNodes,Err)
  !Finish the decomposition
  CALL CMISSDecomposition_CreateFinish(InterfaceDecomposition,Err)

  !Start to create a default (geometric) field on the Interface
  PRINT *, ' == >> CREATING INTERFACE GEOMETRIC FIELD << == '
  CALL CMISSField_Initialise(InterfaceGeometricField,Err)
  CALL CMISSField_CreateStart(InterfaceGeometricFieldUserNumber,Interface,InterfaceGeometricField,Err)
  !Set the decomposition to use
  CALL CMISSField_MeshDecompositionSet(InterfaceGeometricField,InterfaceDecomposition,Err)
  !Set the domain to be used by the field components.
  CALL CMISSField_ComponentMeshComponentSet(InterfaceGeometricField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,Err)
  CALL CMISSField_ComponentMeshComponentSet(InterfaceGeometricField,CMISS_FIELD_U_VARIABLE_TYPE,2,1,Err)
  CALL CMISSField_ComponentMeshComponentSet(InterfaceGeometricField,CMISS_FIELD_U_VARIABLE_TYPE,3,1,Err)
  !Finish creating the first field
  CALL CMISSField_CreateFinish(InterfaceGeometricField,Err)

  !Update the geometric field parameters for the interface field
  PRINT *, ' == >> INITIALIZE INTERFACE GEOMETRIC FIELD FROM INTERFACE GENERATED MESH << == '
  CALL CMISSGeneratedMesh_GeometricParametersCalculate(InterfaceGeneratedMesh,InterfaceGeometricField,Err)

  !Create an interface condition between the two meshes
  PRINT *, ' == >> CREATING INTERFACE CONDITIONS << == '
  CALL CMISSInterfaceCondition_Initialise(InterfaceCondition,Err)
  CALL CMISSInterfaceCondition_CreateStart(InterfaceConditionUserNumber,Interface,InterfaceGeometricField,InterfaceCondition,Err)
  !Specify the method for the interface condition
  IF(PenaltyMethod) THEN
    CALL CMISSInterfaceCondition_MethodSet(InterfaceCondition,CMISS_INTERFACE_CONDITION_PENALTY_METHOD,Err)
  ELSE
    CALL CMISSInterfaceCondition_MethodSet(InterfaceCondition,CMISS_INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,Err)
  ENDIF

  !Specify the type of interface condition operator
  CALL CMISSInterfaceCondition_OperatorSet(InterfaceCondition,CMISS_INTERFACE_CONDITION_FIELD_CONTINUITY_OPERATOR,Err)
  !Add in the dependent variables from the equations sets
  CALL CMISSInterfaceCondition_DependentVariableAdd(InterfaceCondition,Mesh1Index,EquationsSet1, &
    & CMISS_FIELD_U_VARIABLE_TYPE,Err)
  CALL CMISSInterfaceCondition_DependentVariableAdd(InterfaceCondition,Mesh2Index,EquationsSet2, &
    & CMISS_FIELD_U_VARIABLE_TYPE,Err)
  !Finish creating the interface condition
  CALL CMISSInterfaceCondition_CreateFinish(InterfaceCondition,Err)

  !Create the Lagrange multipliers field
  PRINT *, ' == >> CREATING INTERFACE LAGRANGE FIELD << == '
  CALL CMISSField_Initialise(LagrangeField,Err)
  CALL CMISSInterfaceCondition_LagrangeFieldCreateStart(InterfaceCondition,LagrangeFieldUserNumber,LagrangeField,Err)
  !Finish the Lagrange multipliers field
  CALL CMISSInterfaceCondition_LagrangeFieldCreateFinish(InterfaceCondition,Err)

  IF(PenaltyMethod) THEN
    !Create the PenaltyMethod field
    PRINT *, ' == >> CREATING INTERFACE PenaltyMethod FIELD << == '
    CALL CMISSField_Initialise(PenaltyField,Err)
    CALL CMISSInterfaceCondition_PenaltyFieldCreateStart(InterfaceCondition,PenaltyFieldUserNumber, &
      & PenaltyField,Err)
    !Finish the PenaltyMethod field
    CALL CMISSInterfaceCondition_PenaltyFieldCreateFinish(InterfaceCondition,Err)
    !Set the PenaltyMethod field coefficients
    CALL CMISSField_ComponentValuesInitialise(PenaltyField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1, &
      & 100.0_CMISSDP,Err)
    CALL CMISSField_ComponentValuesInitialise(PenaltyField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2, &
      & 100.0_CMISSDP,Err)
    CALL CMISSField_ComponentValuesInitialise(PenaltyField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,3, &
      & 100.0_CMISSDP,Err)
  ENDIF

  !Create the interface condition equations
  PRINT *, ' == >> CREATING INTERFACE EQUATIONS << == '
  CALL CMISSInterfaceEquations_Initialise(InterfaceEquations,Err)
  CALL CMISSInterfaceCondition_EquationsCreateStart(InterfaceCondition,InterfaceEquations,Err)
  !Set the interface equations sparsity
  !CALL CMISSInterfaceEquations_SparsitySet(InterfaceEquations,CMISS_EQUATIONS_SPARSE_MATRICES,Err)
  CALL CMISSInterfaceEquations_SparsitySet(InterfaceEquations,CMISS_EQUATIONS_FULL_MATRICES,Err)
  !Set the interface equations output
  CALL CMISSInterfaceEquations_OutputTypeSet(InterfaceEquations,CMISS_EQUATIONS_NO_OUTPUT,Err)
  !CALL CMISSInterfaceEquations_OutputTypeSet(InterfaceEquations,CMISS_EQUATIONS_TIMING_OUTPUT,Err)
  !CALL CMISSInterfaceEquations_OutputTypeSet(InterfaceEquations,CMISS_EQUATIONS_MATRIX_OUTPUT,Err)
  !CALL CMISSInterfaceEquations_OutputTypeSet(InterfaceEquations,CMISS_EQUATIONS_ELEMENT_MATRIX_OUTPUT,Err)

  !Finish creating the interface equations
  CALL CMISSInterfaceCondition_EquationsCreateFinish(InterfaceCondition,Err)

  !######################################### FINISH CREATING INTERFACE ##########################################

  !Start the creation of a coupled problem.
  PRINT *, ' == >> CREATING PROBLEM << == '
  CALL CMISSProblem_Initialise(CoupledProblem,Err)
  CALL CMISSProblem_CreateStart(CoupledProblemUserNumber,CoupledProblem,Err)
  CALL CMISSProblem_SpecificationSet(CoupledProblem,CMISS_PROBLEM_ELASTICITY_CLASS, &
    & CMISS_PROBLEM_FINITE_ELASTICITY_TYPE,CMISS_PROBLEM_NO_SUBTYPE,Err)
  !Finish the creation of a problem.
  CALL CMISSProblem_CreateFinish(CoupledProblem,Err)

  !Start the creation of the problem control loop for the coupled problem
  PRINT *, ' == >> CREATING PROBLEM CONTROL LOOP << == '
  CALL CMISSProblem_ControlLoopCreateStart(CoupledProblem,Err)
  CALL CMISSControlLoop_Initialise(ControlLoop,Err)
  CALL CMISSProblem_ControlLoopGet(CoupledProblem,CMISS_CONTROL_LOOP_NODE,ControlLoop,Err)
  CALL CMISSControlLoop_MaximumIterationsSet(ControlLoop,LoadSteps,Err)
  !Finish creating the problem control loop
  CALL CMISSProblem_ControlLoopCreateFinish(CoupledProblem,Err)
 
  !Start the creation of the problem solver for the coupled problem
  PRINT *, ' == >> CREATING PROBLEM SOLVERS << == '
  CALL CMISSSolver_Initialise(NonLinearSolver,Err)
  CALL CMISSSolver_Initialise(LinearSolver,Err)
  CALL CMISSProblem_SolversCreateStart(CoupledProblem,Err)
  CALL CMISSProblem_SolverGet(CoupledProblem,CMISS_CONTROL_LOOP_NODE,1,NonLinearSolver,Err)
  !CALL CMISSSolver_OutputTypeSet(NonLinearSolver,CMISS_SOLVER_NO_OUTPUT,Err)
  CALL CMISSSolver_OutputTypeSet(NonLinearSolver,CMISS_SOLVER_PROGRESS_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(NonLinearSolver,CMISS_SOLVER_TIMING_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(NonLinearSolver,CMISS_SOLVER_SOLVER_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(NonLinearSolver,CMISS_SOLVER_MATRIX_OUTPUT,Err)
  CALL CMISSSolver_NewtonJacobianCalculationTypeSet(NonLinearSolver,CMISS_SOLVER_NEWTON_JACOBIAN_FD_CALCULATED,Err)
  !CALL CMISSSolver_NewtonJacobianCalculationTypeSet(NonLinearSolver,CMISS_SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED,Err)
  CALL CMISSSolver_NewtonMaximumFunctionEvaluationsSet(NonLinearSolver,10000,Err)
  !CALL CMISSSolver_NewtonMaximumIterationsSet(NonLinearSolver,5,Err)
  !CALL CMISSSolver_NewtonAbsoluteToleranceSet(NonlinearSolver,1e-6_CMISSDP,Err)
  !CALL CMISSSolver_NewtonRelativeToleranceSet(NonlinearSolver,1e-6_CMISSDP,Err)
  !CALL CMISSSolver_LibraryTypeSet(NonLinearSolver,CMISS_SOLVER_MUMPS_LIBRARY,Err)
  CALL CMISSSolver_NewtonLinearSolverGet(NonLinearSolver,LinearSolver,Err)
  CALL CMISSSolver_LinearTypeSet(LinearSolver,CMISS_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
  CALL CMISSSolver_LibraryTypeSet(LinearSolver,CMISS_SOLVER_LAPACK_LIBRARY,Err)
  !CALL CMISSSolver_LinearTypeSet(LinearSolver,CMISS_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE,Err)
  !CALL CMISSSolver_LinearIterativeMaximumIterationsSet(LinearSolver,1,Err)
  CALL CMISSSolver_OutputTypeSet(LinearSolver,CMISS_SOLVER_PROGRESS_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(LinearSolver,CMISS_SOLVER_TIMING_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(LinearSolver,CMISS_SOLVER_SOLVER_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(LinearSolver,CMISS_SOLVER_MATRIX_OUTPUT,Err)
  !Finish the creation of the problem solver
  CALL CMISSProblem_SolversCreateFinish(CoupledProblem,Err)

  !Start the creation of the problem solver equations for the coupled problem
  PRINT *, ' == >> CREATING PROBLEM SOLVER EQUATIONS << == '
  CALL CMISSSolver_Initialise(NonLinearSolver,Err)
  CALL CMISSSolverEquations_Initialise(NonLinearSolverEquations,Err)
  CALL CMISSProblem_SolverEquationsCreateStart(CoupledProblem,Err)
  !Get the solve equations
  CALL CMISSProblem_SolverGet(CoupledProblem,CMISS_CONTROL_LOOP_NODE,1,NonLinearSolver,Err)
  CALL CMISSSolver_SolverEquationsGet(NonLinearSolver,NonLinearSolverEquations,Err)
  !Set the solver equations sparsity
  !CALL CMISSSolverEquations_SparsityTypeSet(NonLinearSolverEquations,CMISS_SOLVER_SPARSE_MATRICES,Err)
  CALL CMISSSolverEquations_SparsityTypeSet(NonLinearSolverEquations,CMISS_SOLVER_FULL_MATRICES,Err)  
  !Add in the first equations set
  CALL CMISSSolverEquations_EquationsSetAdd(NonLinearSolverEquations,EquationsSet1,EquationsSet1Index,Err)
  !Add in the second equations set
  CALL CMISSSolverEquations_EquationsSetAdd(NonLinearSolverEquations,EquationsSet2,EquationsSet2Index,Err)
  IF(Uncoupled.EQV..FALSE.) THEN
    !Add in the interface condition
    CALL CMISSSolverEquations_InterfaceConditionAdd(NonLinearSolverEquations,InterfaceCondition,InterfaceConditionIndex,Err)
  ENDIF
  !Finish the creation of the problem solver equations
  CALL CMISSProblem_SolverEquationsCreateFinish(CoupledProblem,Err)

  !Export the fields
  PRINT *, ' == >> EXPORTING FIELDS << == '
  CALL CMISSFields_Initialise(Fields1,Err)
  CALL CMISSFields_Create(Region1,Fields1,Err)
  CALL CMISSFields_NodesExport(Fields1,"3DCoupled-FiniteElasticity-FiniteElasticity_1_Initial","FORTRAN",Err)
  CALL CMISSFields_ElementsExport(Fields1,"3DCoupled-FiniteElasticity-FiniteElasticity_1_Initial","FORTRAN",Err)
  CALL CMISSFields_Finalise(Fields1,Err)
  CALL CMISSFields_Initialise(Fields2,Err)
  CALL CMISSFields_Create(Region2,Fields2,Err)
  CALL CMISSFields_NodesExport(Fields2,"3DCoupled-FiniteElasticity-FiniteElasticity_2_Initial","FORTRAN",Err)
  CALL CMISSFields_ElementsExport(Fields2,"3DCoupled-FiniteElasticity-FiniteElasticity_2_Initial","FORTRAN",Err)
  CALL CMISSFields_Finalise(Fields2,Err)
  IF(Uncoupled.EQV..FALSE.) THEN
    CALL CMISSFields_Initialise(Fields3,Err)
    CALL CMISSFields_Create(Interface,Fields3,Err)
    CALL CMISSFields_NodesExport(Fields3,"3DCoupled-FiniteElasticity-FiniteElasticity_Interface_Initial","FORTRAN",Err)
    CALL CMISSFields_ElementsExport(Fields3,"3DCoupled-FiniteElasticity-FiniteElasticity_Interface_Initial","FORTRAN",Err)
    CALL CMISSFields_Finalise(Fields3,Err)
  ENDIF

  !Prescribe boundary conditions (absolute nodal parameters)
  CALL CMISSBoundaryConditions_Initialise(BoundaryConditions,Err)
  CALL CMISSSolverEquations_BoundaryConditionsCreateStart(NonLinearSolverEquations,BoundaryConditions,Err)

  IF(Uncoupled.EQV..FALSE.) THEN
    !============================================================================================================================
    ! COUPLED
    !============================================================================================================================
    !Start the creation of the solver equations boundary conditions for the first region
    PRINT *, ' == >> SETTING BOUNDARY CONDITIONS(1) << == '

    BackNodeGroup = [1,3,5,7]
    DO node_idx=1,SIZE(BackNodeGroup)
      node=BackNodeGroup(node_idx)
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField1,CMISS_FIELD_U_VARIABLE_TYPE,1, &
        & CMISS_NO_GLOBAL_DERIV,node,1,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
      IF(DisplacementInterpolationType==CMISS_BASIS_CUBIC_HERMITE_INTERPOLATION)THEN
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField1,CMISS_FIELD_U_VARIABLE_TYPE,1, &
          & CMISS_GLOBAL_DERIV_S2,node,1,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField1,CMISS_FIELD_U_VARIABLE_TYPE,1, &
          & CMISS_GLOBAL_DERIV_S3,node,1,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField1,CMISS_FIELD_U_VARIABLE_TYPE,1, &
          & CMISS_GLOBAL_DERIV_S2_S3,node,1,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
      ENDIF
    ENDDO

    Side4NodeGroup = [1,2,5,6]
    DO node_idx=1,SIZE(Side4NodeGroup)
      node=Side4NodeGroup(node_idx)
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField1,CMISS_FIELD_U_VARIABLE_TYPE,1, &
        & CMISS_NO_GLOBAL_DERIV,node,2,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
    ENDDO

    Bottom4NodeGroup = [1,2,3,4]
    DO node_idx=1,SIZE(Bottom4NodeGroup)
      node=Bottom4NodeGroup(node_idx)
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField1,CMISS_FIELD_U_VARIABLE_TYPE,1, &
        & CMISS_NO_GLOBAL_DERIV,node,3,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
    ENDDO

    Side4NodeGroup = [1,2,5,6]
    DO node_idx=1,SIZE(Side4NodeGroup)
      node=Side4NodeGroup(node_idx)
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField2,CMISS_FIELD_U_VARIABLE_TYPE,1, &
        & CMISS_NO_GLOBAL_DERIV,node,2,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
    ENDDO

    Bottom4NodeGroup = [1,2,3,4]
    DO node_idx=1,SIZE(Bottom4NodeGroup)
      node=Bottom4NodeGroup(node_idx)
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField2,CMISS_FIELD_U_VARIABLE_TYPE,1, &
        & CMISS_NO_GLOBAL_DERIV,node,3,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
    ENDDO

    FrontNodeGroup = [2,4,6,8]
    DO node_idx=1,SIZE(FrontNodeGroup)
      node=FrontNodeGroup(node_idx)
      IF(ForceBoundaryConditions.EQV..TRUE.) THEN
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField2,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,1, &
          & CMISS_NO_GLOBAL_DERIV,node,1,CMISS_BOUNDARY_CONDITION_FIXED,-1.5_CMISSDP,Err)
      ELSE
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField2,CMISS_FIELD_U_VARIABLE_TYPE,1, &
          & CMISS_NO_GLOBAL_DERIV,node,1,CMISS_BOUNDARY_CONDITION_FIXED,2.1_CMISSDP,Err)
      ENDIF
      IF(DisplacementInterpolationType==CMISS_BASIS_CUBIC_HERMITE_INTERPOLATION)THEN
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField2,CMISS_FIELD_U_VARIABLE_TYPE,1, &
          & CMISS_GLOBAL_DERIV_S2,node,1,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField2,CMISS_FIELD_U_VARIABLE_TYPE,1, &
          & CMISS_GLOBAL_DERIV_S3,node,1,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField2,CMISS_FIELD_U_VARIABLE_TYPE,1, &
          & CMISS_GLOBAL_DERIV_S2_S3,node,1,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
      ENDIF
    ENDDO

    PRINT *, ' == >> SETTING INTERFACE BOUNDARY CONDITIONS << == '
    CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,LagrangeField,CMISS_FIELD_U_VARIABLE_TYPE,1, &
      & CMISS_NO_GLOBAL_DERIV,1,2,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
    CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,LagrangeField,CMISS_FIELD_U_VARIABLE_TYPE,1, &
      & CMISS_NO_GLOBAL_DERIV,3,2,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
    CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,LagrangeField,CMISS_FIELD_U_VARIABLE_TYPE,1, &
      & CMISS_NO_GLOBAL_DERIV,1,3,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
    CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,LagrangeField,CMISS_FIELD_U_VARIABLE_TYPE,1, &
      & CMISS_NO_GLOBAL_DERIV,2,3,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)

    DO node=1,4 !Number of nodes in interface mesh
      IF(DisplacementInterpolationType==CMISS_BASIS_CUBIC_HERMITE_INTERPOLATION)THEN
        DO component=1,3
          CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,LagrangeField,CMISS_FIELD_U_VARIABLE_TYPE,1, &
            & CMISS_GLOBAL_DERIV_S1,node,component,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
          CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,LagrangeField,CMISS_FIELD_U_VARIABLE_TYPE,1, &
            & CMISS_GLOBAL_DERIV_S2,node,component,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
          CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,LagrangeField,CMISS_FIELD_U_VARIABLE_TYPE,1, &
            & CMISS_GLOBAL_DERIV_S1_S2,node,component,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
        ENDDO
      ENDIF
      IF(HydrostaticPressureBasis) THEN
        IF ((Incompressible1.eqv..TRUE.) .AND. (Incompressible2.eqv..TRUE.)) THEN
          component=4
          CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,LagrangeField,CMISS_FIELD_U_VARIABLE_TYPE,1, &
            & CMISS_NO_GLOBAL_DERIV,node,component,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
          IF(DisplacementInterpolationType==CMISS_BASIS_CUBIC_HERMITE_INTERPOLATION)THEN
            CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,LagrangeField,CMISS_FIELD_U_VARIABLE_TYPE,1, &
              & CMISS_GLOBAL_DERIV_S1,node,component,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
            CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,LagrangeField,CMISS_FIELD_U_VARIABLE_TYPE,1, &
              & CMISS_GLOBAL_DERIV_S2,node,component,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
            CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,LagrangeField,CMISS_FIELD_U_VARIABLE_TYPE,1, &
              & CMISS_GLOBAL_DERIV_S1_S2,node,component,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
          ENDIF
        ENDIF
      ENDIF 
    ENDDO

  ELSE
    !============================================================================================================================
    ! Uncoupled
    !============================================================================================================================
    !Start the creation of the solver equations boundary conditions for the first region
    PRINT *, ' == >> SETTING BOUNDARY CONDITIONS(1) << == '

    BackNodeGroup = [1,4,7,10]
    DO node_idx=1,SIZE(BackNodeGroup)
      node=BackNodeGroup(node_idx)
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField1,CMISS_FIELD_U_VARIABLE_TYPE,1, &
        & CMISS_NO_GLOBAL_DERIV,node,1,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
      IF(DisplacementInterpolationType==CMISS_BASIS_CUBIC_HERMITE_INTERPOLATION)THEN
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField1,CMISS_FIELD_U_VARIABLE_TYPE,1, &
          & CMISS_GLOBAL_DERIV_S2,node,1,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField1,CMISS_FIELD_U_VARIABLE_TYPE,1, &
          & CMISS_GLOBAL_DERIV_S3,node,1,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField1,CMISS_FIELD_U_VARIABLE_TYPE,1, &
          & CMISS_GLOBAL_DERIV_S2_S3,node,1,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
      ENDIF
    ENDDO

    FrontNodeGroup = [3,6,9,12]
    DO node_idx=1,SIZE(FrontNodeGroup)
      node=FrontNodeGroup(node_idx)
      IF(ForceBoundaryConditions.EQV..TRUE.) THEN
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField1,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,1, &
          & CMISS_NO_GLOBAL_DERIV,node,1,CMISS_BOUNDARY_CONDITION_FIXED,-1.5_CMISSDP,Err)
      ELSE
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField1,CMISS_FIELD_U_VARIABLE_TYPE,1, &
          & CMISS_NO_GLOBAL_DERIV,node,1,CMISS_BOUNDARY_CONDITION_FIXED,2.1_CMISSDP,Err)
      ENDIF
      IF(DisplacementInterpolationType==CMISS_BASIS_CUBIC_HERMITE_INTERPOLATION)THEN
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField1,CMISS_FIELD_U_VARIABLE_TYPE,1, &
          & CMISS_GLOBAL_DERIV_S2,node,1,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField1,CMISS_FIELD_U_VARIABLE_TYPE,1, &
          & CMISS_GLOBAL_DERIV_S3,node,1,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField1,CMISS_FIELD_U_VARIABLE_TYPE,1, &
          & CMISS_GLOBAL_DERIV_S2_S3,node,1,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
      ENDIF
    ENDDO

    Side6NodeGroup = [1,2,3,7,8,9]
    DO node_idx=1,SIZE(Side6NodeGroup)
      node=Side6NodeGroup(node_idx)
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField1,CMISS_FIELD_U_VARIABLE_TYPE,1, &
        & CMISS_NO_GLOBAL_DERIV,node,2,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
    ENDDO

    Bottom6NodeGroup = [1,2,3,4,5,6]
    DO node_idx=1,SIZE(Bottom6NodeGroup)
      node=Bottom6NodeGroup(node_idx)
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField1,CMISS_FIELD_U_VARIABLE_TYPE,1, &
        & CMISS_NO_GLOBAL_DERIV,node,3,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
    ENDDO

    !Start the creation of the solver equations boundary conditions for the second region
    PRINT *, ' == >> SETTING BOUNDARY CONDITIONS(2) << == '

    BackNodeGroup = [1,4,7,10]
    DO node_idx=1,SIZE(BackNodeGroup)
      node=BackNodeGroup(node_idx)
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField2,CMISS_FIELD_U_VARIABLE_TYPE,1, &
        & CMISS_NO_GLOBAL_DERIV,node,1,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
      IF(DisplacementInterpolationType==CMISS_BASIS_CUBIC_HERMITE_INTERPOLATION)THEN
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField2,CMISS_FIELD_U_VARIABLE_TYPE,1, &
          & CMISS_GLOBAL_DERIV_S2,node,1,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField2,CMISS_FIELD_U_VARIABLE_TYPE,1, &
          & CMISS_GLOBAL_DERIV_S3,node,1,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField2,CMISS_FIELD_U_VARIABLE_TYPE,1, &
          & CMISS_GLOBAL_DERIV_S2_S3,node,1,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
      ENDIF
    ENDDO

    FrontNodeGroup = [3,6,9,12]
    DO node_idx=1,SIZE(FrontNodeGroup)
      node=FrontNodeGroup(node_idx)
      IF(ForceBoundaryConditions.EQV..TRUE.) THEN
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField2,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,1, &
          & CMISS_NO_GLOBAL_DERIV,node,1,CMISS_BOUNDARY_CONDITION_FIXED,-1.5_CMISSDP,Err)
      ELSE
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField2,CMISS_FIELD_U_VARIABLE_TYPE,1, &
          & CMISS_NO_GLOBAL_DERIV,node,1,CMISS_BOUNDARY_CONDITION_FIXED,2.1_CMISSDP,Err)
      ENDIF
      IF(DisplacementInterpolationType==CMISS_BASIS_CUBIC_HERMITE_INTERPOLATION)THEN
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField2,CMISS_FIELD_U_VARIABLE_TYPE,1, &
          & CMISS_GLOBAL_DERIV_S2,node,1,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField2,CMISS_FIELD_U_VARIABLE_TYPE,1, &
          & CMISS_GLOBAL_DERIV_S3,node,1,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField2,CMISS_FIELD_U_VARIABLE_TYPE,1, &
          & CMISS_GLOBAL_DERIV_S2_S3,node,1,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
      ENDIF
    ENDDO

    Side6NodeGroup = [1,2,3,7,8,9]
    DO node_idx=1,SIZE(Side6NodeGroup)
      node=Side6NodeGroup(node_idx)
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField2,CMISS_FIELD_U_VARIABLE_TYPE,1, &
        & CMISS_NO_GLOBAL_DERIV,node,2,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
    ENDDO

    Bottom6NodeGroup = [1,2,3,4,5,6]
    DO node_idx=1,SIZE(Bottom6NodeGroup)
      node=Bottom6NodeGroup(node_idx)
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField2,CMISS_FIELD_U_VARIABLE_TYPE,1, &
        & CMISS_NO_GLOBAL_DERIV,node,3,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
    ENDDO

  ENDIF

  CALL CMISSSolverEquations_BoundaryConditionsCreateFinish(NonLinearSolverEquations,Err)

  !Solve the problem
  PRINT *, ' == >> SOLVING PROBLEM << == '
  CALL CMISSProblem_Solve(CoupledProblem,Err)

  !Export the fields
  PRINT *, ' == >> EXPORTING FIELDS << == '
  CALL CMISSFields_Initialise(Fields1,Err)
  CALL CMISSFields_Create(Region1,Fields1,Err)
  CALL CMISSFields_NodesExport(Fields1,"3DCoupled-FiniteElasticity-FiniteElasticity_1","FORTRAN",Err)
  CALL CMISSFields_ElementsExport(Fields1,"3DCoupled-FiniteElasticity-FiniteElasticity_1","FORTRAN",Err)
  CALL CMISSFields_Finalise(Fields1,Err)
  CALL CMISSFields_Initialise(Fields2,Err)
  CALL CMISSFields_Create(Region2,Fields2,Err)
  CALL CMISSFields_NodesExport(Fields2,"3DCoupled-FiniteElasticity-FiniteElasticity_2","FORTRAN",Err)
  CALL CMISSFields_ElementsExport(Fields2,"3DCoupled-FiniteElasticity-FiniteElasticity_2","FORTRAN",Err)
  CALL CMISSFields_Finalise(Fields2,Err)
  IF(Uncoupled.EQV..FALSE.) THEN
    CALL CMISSFields_Initialise(Fields3,Err)
    CALL CMISSFields_Create(Interface,Fields3,Err)
    CALL CMISSFields_NodesExport(Fields3,"3DCoupled-FiniteElasticity-FiniteElasticity_Interface","FORTRAN",Err)
    CALL CMISSFields_ElementsExport(Fields3,"3DCoupled-FiniteElasticity-FiniteElasticity_Interface","FORTRAN",Err)
    CALL CMISSFields_Finalise(Fields3,Err)
  ENDIF
  
  !Finialise CMISS
  CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."

  STOP
  
CONTAINS

  SUBROUTINE HANDLE_ERROR(ERROR_STRING)

    CHARACTER(LEN=*), INTENT(IN) :: ERROR_STRING

    WRITE(*,'(">>ERROR: ",A)') ERROR_STRING(1:LEN_TRIM(ERROR_STRING))
    STOP

  END SUBROUTINE HANDLE_ERROR

END PROGRAM COUPLEDFINITEELASTICITYFINITEELASTICITY
