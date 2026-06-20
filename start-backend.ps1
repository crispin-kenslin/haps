$ErrorActionPreference = 'Stop'

$scriptRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $scriptRoot

$pythonExe = Join-Path $scriptRoot 'herb-web\Scripts\python.exe'
if (-not (Test-Path $pythonExe)) {
    throw 'Python virtual environment not found at herb-web\Scripts\python.exe'
}

$logDir = Join-Path $scriptRoot 'logs'
if (-not (Test-Path $logDir)) {
    New-Item -ItemType Directory -Path $logDir -Force | Out-Null
}

$logFile = Join-Path $logDir 'backend-startup.log'

Add-Content -Path $logFile -Value "[$(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')] Starting app.py"

$ErrorActionPreference = 'Continue'

# Run app.py and log all output
& $pythonExe app.py *>> $logFile

$exitCode = $LASTEXITCODE

Add-Content -Path $logFile -Value "[$(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')] app.py exited with code $exitCode"

if ($exitCode -ne 0) {
    exit $exitCode
}