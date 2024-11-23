# IMPORTANT:
# 1. move to the correct top level directory
# 2. make sure the input files are edited correctly, e.g. pitch angle...

# Get the current directory
$currentDir = Get-Location

# Define the source and target directories
$sourceDir = Join-Path $currentDir "Lab_Files\panel_code\Windows_Executables"
$targetBaseDir = Join-Path $currentDir "discretisation"

# Ensure the target base directory exists
if (!(Test-Path $targetBaseDir)) {
    New-Item -ItemType Directory -Path $targetBaseDir
}

# Path to the log file
$logFile = "runtimes.log"

# Ensure the log file is empty or create it if it doesn't exist
Set-Content -Path $logFile -Value ""

# Loop through m (1 to 30) and n (1 to 20)
for ($m = 1; $m -le 30; $m++) {
    for ($n = 1; $n -le 20; $n++) {
        # Construct folder name c$(m)-s$(n)
        $folderName = "c$($m)-s$($n)"
        $folderPath = Join-Path $targetBaseDir $folderName

        # Create the folder if it doesn't exist
        if (!(Test-Path $folderPath)) {
            New-Item -ItemType Directory -Path $folderPath
        }

        # Copy files from source to the newly created folder
        Copy-Item -Path "$sourceDir\*" -Destination $folderPath -Recurse -Force

        cd $folderPath
        # Define the path to the Fortran executable
        $fortranExecutable = Join-Path -Path $folderPath -ChildPath "geowing.exe"

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

        # measure the computation time for running the panel code
        $runtime = Measure-Command {

        # Define the path to the Fortran executable
        $fortranExecutable2 = Join-Path -Path $folderPath -ChildPath "combo.exe"

        # Run the Fortran executable with the input data piped to it
        & $fortranExecutable2

        }

        $runtime = $runtime.TotalSeconds

        # output the computation time
        Add-Content -Path $logFile -Value $runtime

        Write-Host "Fortran executable 2 has finished running."






    }
}

Write-Host "Folders and file copies completed successfully!"
