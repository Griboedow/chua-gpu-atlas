param(
    [string]$ToolsRoot = 'D:\_TMP\chua-proof-tools'
)
$ErrorActionPreference = 'Stop'
$repository = Split-Path -Parent $PSScriptRoot
$tools = [IO.Path]::GetFullPath($ToolsRoot)
if ([IO.Path]::GetPathRoot($repository) -ne 'D:\' -or [IO.Path]::GetPathRoot($tools) -ne 'D:\') {
    throw 'Run this script from the research copy on D: and keep the proof tools on D:.'
}
$python = Join-Path $repository '.venv\Scripts\python.exe'
$elan = Join-Path $tools 'elan'
$lake = Join-Path $elan 'bin\lake.exe'
if (!(Test-Path -LiteralPath $python) -or !(Test-Path -LiteralPath $lake)) {
    throw 'The local Python environment or Lean installation is missing; see formal/README.md.'
}
$temporary = Join-Path $repository '.publication\local-proof-temp'
New-Item -ItemType Directory -Path $temporary -Force | Out-Null
$settings = @{
    ELAN_HOME = $elan
    PATH = (Join-Path $elan 'bin') + [IO.Path]::PathSeparator + $env:PATH
    MATHLIB_CACHE_DIR = (Join-Path $tools 'mathlib-cache')
    PIP_CACHE_DIR = 'D:\_TMP\python-cache\pip'
    TEMP = $temporary
    TMP = $temporary
    PYTHONDONTWRITEBYTECODE = '1'
}
$previous = @{}
foreach ($name in $settings.Keys) {
    $previous[$name] = [Environment]::GetEnvironmentVariable($name, 'Process')
    [Environment]::SetEnvironmentVariable($name, $settings[$name], 'Process')
}
Push-Location -LiteralPath $repository
try {
    & $python -X utf8 -B research/check_chua_formal.py
    if ($LASTEXITCODE -ne 0) { throw "Local proof check failed with exit code $LASTEXITCODE." }
} finally {
    Pop-Location
    foreach ($name in $previous.Keys) {
        [Environment]::SetEnvironmentVariable($name, $previous[$name], 'Process')
    }
}
