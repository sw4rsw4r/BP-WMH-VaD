$ErrorActionPreference = 'Stop'
$Root = Split-Path -Parent $PSScriptRoot
$KgDir = Join-Path $Root 'data\raw\1kg'
$PlinkDir = Join-Path $Root 'software\plink2'
$Log = Join-Path $Root 'logs\download_ld_plink.log'

function Get-PublicFile {
    param([string]$Url, [string]$Destination)
    if (Test-Path -LiteralPath $Destination) {
        Write-Output "SKIP existing: $Destination"
        return
    }
    $partial = "$Destination.part"
    $args = @('-L', '--fail', '--retry', '8', '--retry-delay', '5', '--silent', '--show-error')
    if (Test-Path -LiteralPath $partial) { $args += @('--continue-at', '-') }
    $args += @('--output', $partial, $Url)
    Write-Output "DOWNLOAD: $Url"
    & curl.exe @args
    if ($LASTEXITCODE -ne 0) { throw "curl failed ($LASTEXITCODE): $Url" }
    Move-Item -LiteralPath $partial -Destination $Destination -Force
    Write-Output "DONE: $Destination"
}

Start-Transcript -Path $Log -Append | Out-Null
try {
    New-Item -ItemType Directory -Path $KgDir -Force | Out-Null
    New-Item -ItemType Directory -Path $PlinkDir -Force | Out-Null
    $downloads = @(
        [pscustomobject]@{Url='https://www.dropbox.com/s/y6ytfoybz48dc0u/all_phase3.pgen.zst?dl=1'; Out=(Join-Path $KgDir 'all_phase3.pgen.zst')},
        [pscustomobject]@{Url='https://www.dropbox.com/s/c95n8quqwqww4s0/all_phase3_noannot.pvar.zst?dl=1'; Out=(Join-Path $KgDir 'all_phase3.pvar.zst')},
        [pscustomobject]@{Url='https://www.dropbox.com/scl/fi/haqvrumpuzfutklstazwk/phase3_corrected.psam?dl=1&rlkey=0yyifzj2fb863ddbmsv4jkeq6'; Out=(Join-Path $KgDir 'all_phase3.psam')},
        [pscustomobject]@{Url='https://s3.amazonaws.com/plink2-assets/alpha7/plink2_win_avx2_20260818.zip'; Out=(Join-Path $PlinkDir 'plink2_win_avx2_20260818.zip')}
    )
    foreach ($item in $downloads) { Get-PublicFile -Url $item.Url -Destination $item.Out }

    $plink = Join-Path $PlinkDir 'plink2.exe'
    if (-not (Test-Path -LiteralPath $plink)) {
        Expand-Archive -LiteralPath (Join-Path $PlinkDir 'plink2_win_avx2_20260818.zip') -DestinationPath $PlinkDir -Force
    }
    if (-not (Test-Path -LiteralPath $plink)) { throw "Missing PLINK executable: $plink" }

    $pgen = Join-Path $KgDir 'all_phase3.pgen'
    if (-not (Test-Path -LiteralPath $pgen)) {
        Write-Output 'DECOMPRESS: all_phase3.pgen.zst'
        & $plink --zst-decompress (Join-Path $KgDir 'all_phase3.pgen.zst') $pgen
        if ($LASTEXITCODE -ne 0) { throw 'PLINK zstd decompression failed' }
    }
    if (-not (Test-Path -LiteralPath $pgen)) { throw "Missing decompressed PGEN: $pgen" }

    $records = foreach ($item in $downloads) {
        $f = Get-Item -LiteralPath $item.Out
        $h = Get-FileHash -LiteralPath $item.Out -Algorithm SHA256
        [pscustomobject]@{file=$f.Name; bytes=$f.Length; sha256=$h.Hash.ToLower(); url=$item.Url}
    }
    $records | Export-Csv -LiteralPath (Join-Path $Root 'data\raw\ld_software_download_manifest.tsv') -Delimiter "`t" -NoTypeInformation
    Write-Output (& $plink --version 2>&1 | Select-Object -First 1)
    Write-Output 'LD reference and PLINK ready.'
} finally {
    Stop-Transcript | Out-Null
}
