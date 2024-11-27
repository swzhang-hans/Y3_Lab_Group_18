# Define the path to the Fortran executable
$fortranExecutable = "C:\Users\yasin\Desktop\Coding Projects\Y3_Lab_Group_18\Incompressible\Lab_Files\panel_code\Windows_Executables\geowing.exe"

# Define the input data for the Fortran executable
# Adjust the values as needed for your case
$aerofoil = "N"
$chordwise = 16
$rootChord = 0.225
$semiSpan = 0.6
$taperRatio = 0.42
$sweep = "sweep"
$leadingEdge = 36.9
$dihedral = 0
$twist = -1.5
$rootInc = 0
$spanwise = 9

# Prepare the input data as a string
$inputData = @"
$aerofoil
$chordwise
$rootChord
$semiSpan
$taperRatio
$sweep
$leadingEdge
$dihedral
$twist
$rootInc
$spanwise
"@

# Run the Fortran executable with the input data piped to it
$inputData | & $fortranExecutable

Write-Host "Fortran executable has finished running."

# Define the path to the Fortran executable
$fortranExecutable2 = "C:\Users\yasin\Desktop\Coding Projects\Y3_Lab_Group_18\Incompressible\Lab_Files\panel_code\Windows_Executables\combo.exe"

# Run the Fortran executable with the input data piped to it
& $fortranExecutable2

Write-Host "Fortran executable 2 has finished running."

