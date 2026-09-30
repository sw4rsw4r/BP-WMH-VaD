$ErrorActionPreference = 'Stop'
$Root = Split-Path -Parent $PSScriptRoot
$BbjDir = Join-Path $Root 'data\raw\bbj'
$JpscDir = Join-Path $Root 'data\raw\jpsc'
$Log = Join-Path $Root 'logs\download_gwas.log'

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
    New-Item -ItemType Directory -Path $BbjDir -Force | Out-Null
    New-Item -ItemType Directory -Path $JpscDir -Force | Out-Null
    $downloads = @(
        [pscustomobject]@{Url='https://humandbs.dbcls.jp/files/hum0014/hum0014.v8.SBP.zip'; Out=(Join-Path $BbjDir 'hum0014.v8.SBP.zip')},
        [pscustomobject]@{Url='https://humandbs.dbcls.jp/files/hum0014/hum0014.v8.DBP.zip'; Out=(Join-Path $BbjDir 'hum0014.v8.DBP.zip')},
        [pscustomobject]@{Url='https://humandbs.dbcls.jp/files/hum0014/hum0014.v8.PP.zip'; Out=(Join-Path $BbjDir 'hum0014.v8.PP.zip')},
        [pscustomobject]@{Url='https://humandbs.dbcls.jp/files/hum0014/hum0014.v8.eGFR.zip'; Out=(Join-Path $BbjDir 'hum0014.v8.eGFR.zip')},
        [pscustomobject]@{Url='https://humandbs.dbcls.jp/files/hum0014/hum0014.v8.sCr.zip'; Out=(Join-Path $BbjDir 'hum0014.v8.sCr.zip')},
        [pscustomobject]@{Url='https://humandbs.dbcls.jp/files/hum0014/hum0014.v8.HbA1c.zip'; Out=(Join-Path $BbjDir 'hum0014.v8.HbA1c.zip')},
        [pscustomobject]@{Url='https://humandbs.dbcls.jp/files/hum0014/hum0014.v8.HDL.zip'; Out=(Join-Path $BbjDir 'hum0014.v8.HDL.zip')},
        [pscustomobject]@{Url='https://humandbs.dbcls.jp/files/hum0014/hum0014.v8.LDL.zip'; Out=(Join-Path $BbjDir 'hum0014.v8.LDL.zip')},
        [pscustomobject]@{Url='https://humandbs.dbcls.jp/files/hum0014/hum0014.v8.TG.zip'; Out=(Join-Path $BbjDir 'hum0014.v8.TG.zip')},
        [pscustomobject]@{Url='https://humandbs.dbcls.jp/files/hum0014/hum0014.v8.CRP.zip'; Out=(Join-Path $BbjDir 'hum0014.v8.CRP.zip')},
        [pscustomobject]@{Url='https://humandbs.dbcls.jp/files/hum0466/hum0466.v1.gwas.v1.zip'; Out=(Join-Path $JpscDir 'hum0466.v1.gwas.v1.zip')}
    )
    foreach ($item in $downloads) { Get-PublicFile -Url $item.Url -Destination $item.Out }

    foreach ($trait in @('SBP','DBP','PP','eGFR','sCr','HbA1c','HDL','LDL','TG','CRP')) {
        $zip = Join-Path $BbjDir "hum0014.v8.$trait.zip"
        $target = Join-Path $BbjDir "hum0014.v7.$trait\BBJ.$trait.autosome.txt"
        if (-not (Test-Path -LiteralPath $target)) {
            Write-Output "EXTRACT: $zip"
            Expand-Archive -LiteralPath $zip -DestinationPath $BbjDir -Force
        }
        if (-not (Test-Path -LiteralPath $target)) { throw "Missing extracted file: $target" }
    }
    $outcome = Join-Path $JpscDir 'jpsc_gwasestimates_wmh_MAF0005_Rsq070_ndbc_v2.txt'
    if (-not (Test-Path -LiteralPath $outcome)) {
        Write-Output 'EXTRACT: JPSC-AD WMH'
        Expand-Archive -LiteralPath (Join-Path $JpscDir 'hum0466.v1.gwas.v1.zip') -DestinationPath $JpscDir -Force
    }
    if (-not (Test-Path -LiteralPath $outcome)) { throw "Missing extracted file: $outcome" }

    $records = foreach ($item in $downloads) {
        $f = Get-Item -LiteralPath $item.Out
        $h = Get-FileHash -LiteralPath $item.Out -Algorithm SHA256
        [pscustomobject]@{file=$f.Name; bytes=$f.Length; sha256=$h.Hash.ToLower(); url=$item.Url}
    }
    $records | Export-Csv -LiteralPath (Join-Path $Root 'data\raw\gwas_download_manifest.tsv') -Delimiter "`t" -NoTypeInformation
    Write-Output 'GWAS inputs ready.'
} finally {
    Stop-Transcript | Out-Null
}
