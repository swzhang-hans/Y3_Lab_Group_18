# Get the current directory
$currentDir = Get-Location

# Define the source and target directories
$sourceDir = Join-Path $currentDir "panel_code\Windows_Executables"
$targetBaseDir = Join-Path $currentDir "discretisation"

# Ensure the target base directory exists
if (!(Test-Path $targetBaseDir)) {
    New-Item -ItemType Directory -Path $targetBaseDir
}

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
    }
}

Write-Host "Folders and file copies completed successfully!"
