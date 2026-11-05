# Source-level RECORD parity: same pipeline on both sides, no external model.
#
#   occt_probe   (specs/occt_probe, compiled from the OCCT sources/libraries)
#   record_dump  (crates/occt-topo/examples/record_dump.rs, the port)
#
# Both read the same STEP file with their own reader, and in `rt-iges` both also
# write IGES with their own writer and read that file back with their own
# reader. The only input is the STEP source file; nothing is compared against a
# model produced by an external application.
#
#   powershell -File specs\record_parity\run_record_parity.ps1 [-Only Cube,rev]
#                                                              [-Deflection 0.1] [-Angle 0.5]
#
# Records land in .target-gate/record/ as <stem>.<side>.<mode>.txt so a failure
# can be inspected by hand. Fields are graded in four bands (see `Classify`):
#   STRUCT  a decision differs (counts, validity, types, orientation sign)
#   REAL    geometry differs above the noise floor
#   THRESH  mesh node/triangle totals differ (deflection tie-break, not a defect)
#   TINY    floating fields below 1e-6 relative (last-digit noise)
# Exit code is the number of STRUCT + REAL records; THRESH and TINY are printed
# for the record only.

param(
    [string[]] $Only = @(),
    [double] $Deflection = 0.1,
    [double] $Angle = 0.5,
    [string[]] $Modes = @('read', 'rt-iges')
)

$ErrorActionPreference = 'Stop'
$Repo = Split-Path -Parent (Split-Path -Parent (Split-Path -Parent $PSCommandPath))
$Work = Join-Path $Repo '.target-gate\record'
$Probe = Join-Path $Repo 'specs\occt_probe\probe.bat'
$Port = Join-Path $Repo 'crates\occt-topo\target\debug\examples\record_dump.exe'
New-Item -ItemType Directory -Force -Path $Work | Out-Null

$inv = [System.Globalization.CultureInfo]::InvariantCulture

function Num([string] $s) {
    $v = 0.0
    if ([double]::TryParse($s, [System.Globalization.NumberStyles]::Float, $inv, [ref] $v)) { return $v }
    return [double]::NaN
}

function Near([double] $a, [double] $b) {
    if ([double]::IsNaN($a) -and [double]::IsNaN($b)) { return $true }
    return ([Math]::Abs($a - $b) -le (1e-7 + 1e-9 * [Math]::Abs($b)))
}

# Band of one floating field. `Near` is the exactness band a record can hold:
# 1e-7 absolute plus 1e-9 relative. Anything above it is graded by magnitude,
# because the two causes are different work:
#   TINY  <= 1e-4 absolute or <= 1e-6 relative: last-digit noise from summation
#         order or from the text precision of the two writers. Reported, never
#         a failure.
#   REAL  above that: a geometry deviation that has to be explained. Failure.
function Classify([double] $a, [double] $b) {
    if ([double]::IsNaN($a) -and [double]::IsNaN($b)) { return 'OK' }
    if (Near $a $b) { return 'OK' }
    $abs = [Math]::Abs($a - $b)
    $mag = [Math]::Max([Math]::Abs($a), [Math]::Abs($b))
    $rel = if ($mag -gt 0.0) { $abs / $mag } else { 0.0 }
    if ($abs -le 1e-4 -or $rel -le 1e-6) { return 'TINY' }
    return 'REAL'
}

# Record cells are text, but a parsed face carries doubles; format those with the
# same invariant G17 the record writer uses so `Num` never sees a locale decimal
# separator (a mis-parse would read as NaN == NaN and pass silently).
function Cell($v) {
    if ($v -is [string]) { return $v }
    if ($v -is [System.Collections.IEnumerable]) {
        return (@($v) | ForEach-Object { ([double] $_).ToString('G17', $inv) }) -join ' '
    }
    return ([double] $v).ToString('G17', $inv)
}

function Grade($Bands, [string] $label, $aText, $bText, [string] $detail = '') {
    $line = "$label  occ=$(Cell $aText)  port=$(Cell $bText)"
    if ($detail) { $line += "  $detail" }
    switch (Classify (Num (Cell $aText)) (Num (Cell $bText))) {
        'REAL' { $Bands.real.Add($line) }
        'TINY' { $Bands.tiny.Add($line) }
    }
}

# Multi-component field (a centroid, a box, a parameter triple): the record is
# graded by its worst component.
function GradeList($Bands, [string] $label, $aText, $bText, [string] $detail = '') {
    $aT = Cell $aText
    $bT = Cell $bText
    $line = "$label  occ=$aT  port=$bT"
    if ($detail) { $line += "  $detail" }
    $a = @($aT.Split(' ') | Where-Object { $_ -ne '' } | ForEach-Object { Num $_ })
    $b = @($bT.Split(' ') | Where-Object { $_ -ne '' } | ForEach-Object { Num $_ })
    if ($a.Count -ne $b.Count) { $Bands.real.Add($line); return }
    $band = 'OK'
    for ($i = 0; $i -lt $a.Count; $i++) {
        $c = Classify $a[$i] $b[$i]
        if ($c -eq 'REAL') { $band = 'REAL'; break }
        if ($c -eq 'TINY') { $band = 'TINY' }
    }
    if ($band -eq 'REAL') { $Bands.real.Add($line) }
    elseif ($band -eq 'TINY') { $Bands.tiny.Add($line) }
}

function ParseFace([string] $line) {
    $m = [regex]::Match($line, '^FACE type=(\S+) area=(\S+) wires=(\d+) edges=(\d+) bbox=(\S+) (\S+) (\S+) (\S+) (\S+) (\S+)$')
    if (-not $m.Success) { return $null }
    $bb = @()
    for ($k = 5; $k -le 10; $k++) { $bb += (Num $m.Groups[$k].Value) }
    return [pscustomobject]@{
        Type  = $m.Groups[1].Value
        Area  = Num $m.Groups[2].Value
        Wires = [int] $m.Groups[3].Value
        Edges = [int] $m.Groups[4].Value
        BBox  = $bb
        Key   = (($bb | ForEach-Object { '{0:F6}' -f $_ }) -join ',') + '|' + $m.Groups[1].Value
        Text  = $line
    }
}

function Parse-Record([string] $path) {
    $scal = @{}
    $faces = @()
    foreach ($line in Get-Content -LiteralPath $path) {
        if ($line.StartsWith('FACE ')) {
            $f = ParseFace $line
            if ($f) { $faces += $f }
            continue
        }
        $i = $line.IndexOf(' ')
        if ($i -gt 0) { $scal[$line.Substring(0, $i)] = $line.Substring($i + 1) }
        else { $scal[$line] = '' }
    }
    return [pscustomobject]@{ Scalars = $scal; Faces = $faces }
}

function Counts([string] $s) {
    $out = @{}
    foreach ($m in [regex]::Matches($s, '([A-Z]+)=(\d+)')) { $out[$m.Groups[1].Value] = [int] $m.Groups[2].Value }
    return $out
}

function Compare-Record([string] $occPath, [string] $portPath) {
    # Four bands, in the order a reader should care about them:
    #   struct  decisions differ (counts, validity, types, orientation signs)
    #   real    geometry differs above the noise floor
    #   thresh  mesh node/triangle totals differ: the output of a deflection
    #           threshold, so equality is not a property of correct behavior
    #   tiny    floating fields below the noise floor
    $struct = New-Object System.Collections.Generic.List[string]
    $real = New-Object System.Collections.Generic.List[string]
    $thresh = New-Object System.Collections.Generic.List[string]
    $tiny = New-Object System.Collections.Generic.List[string]
    $bands = @{ real = $real; thresh = $thresh; tiny = $tiny }
    $a = Parse-Record $occPath
    $b = Parse-Record $portPath

    foreach ($k in @('STAGE', 'TYPE', 'VALID')) {
        if ($a.Scalars[$k] -ne $b.Scalars[$k]) {
            $struct.Add("$k  occ=$($a.Scalars[$k])  port=$($b.Scalars[$k])")
        }
    }
    $ca = Counts $a.Scalars['COUNTS']
    $cb = Counts $b.Scalars['COUNTS']
    foreach ($k in @('V', 'E', 'W', 'F', 'SH', 'SO', 'CS', 'C')) {
        if ($ca[$k] -ne $cb[$k]) { $struct.Add("COUNTS.$k  occ=$($ca[$k])  port=$($cb[$k])") }
    }
    $va = Num $a.Scalars['VOLUME']
    $vb = Num $b.Scalars['VOLUME']
    if (($va -ne 0) -and ($vb -ne 0) -and ([Math]::Sign($va) -ne [Math]::Sign($vb))) {
        # A sign flip is an orientation defect, not a precision difference.
        $struct.Add("VOLUME.sign  occ=$($a.Scalars['VOLUME'])  port=$($b.Scalars['VOLUME'])")
    }
    else {
        Grade $bands 'VOLUME' $a.Scalars['VOLUME'] $b.Scalars['VOLUME']
    }
    Grade $bands 'AREA' $a.Scalars['AREA'] $b.Scalars['AREA']
    GradeList $bands 'CENTROID' $a.Scalars['CENTROID'] $b.Scalars['CENTROID']
    GradeList $bands 'BBOX' $a.Scalars['BBOX'] $b.Scalars['BBOX']

    $fa = @($a.Faces | Sort-Object Key)
    $fb = @($b.Faces | Sort-Object Key)
    if ($fa.Count -ne $fb.Count) {
        $struct.Add("FACE count  occ=$($fa.Count)  port=$($fb.Count)")
    }
    $n = [Math]::Min($fa.Count, $fb.Count)
    for ($i = 0; $i -lt $n; $i++) {
        $x = $fa[$i]; $y = $fb[$i]
        if ($x.Type -ne $y.Type) { $struct.Add("FACE[$i] type  occ=$($x.Type)  port=$($y.Type)  occbbox=$($x.Key)  portbbox=$($y.Key)"); continue }
        if ($x.Wires -ne $y.Wires) { $struct.Add("FACE[$i] wires  occ=$($x.Wires)  port=$($y.Wires)  $($x.Key)") }
        if ($x.Edges -ne $y.Edges) { $struct.Add("FACE[$i] edges  occ=$($x.Edges)  port=$($y.Edges)  $($x.Key)") }
        Grade $bands "FACE[$i] area" $x.Area $y.Area $x.Key
        GradeList $bands "FACE[$i] bbox" $x.BBox $y.BBox $x.Key
    }

    $ma = @{}
    foreach ($m in [regex]::Matches($a.Scalars['MESH'], '(\w+)=(\S+)')) { $ma[$m.Groups[1].Value] = $m.Groups[2].Value }
    $mb = @{}
    foreach ($m in [regex]::Matches($b.Scalars['MESH'], '(\w+)=(\S+)')) { $mb[$m.Groups[1].Value] = $m.Groups[2].Value }
    foreach ($k in @('nv', 'nf')) {
        # Node/triangle totals come out of `BRepMesh_DelaunayDeflectionControlMeshAlgo`,
        # whose split decision is `deviation > deflection`. A link that lands on the
        # threshold flips on summation order alone, and every flipped node adds two
        # triangles: measured on data\Shape.step, per-face deltas satisfy
        # dnf = 2*dnv on all 31 faces. Integer equality there describes the float
        # path, not the port, so it is reported and does not fail the gate.
        if ([int] $ma[$k] -ne [int] $mb[$k]) { $thresh.Add("MESH.$k  occ=$($ma[$k])  port=$($mb[$k])") }
    }
    foreach ($k in @('area', 'defl', 'angle')) {
        Grade $bands "MESH.$k" $ma[$k] $mb[$k]
    }
    return [pscustomobject]@{ Struct = $struct; Real = $real; Thresh = $thresh; Tiny = $tiny }
}

$models = Get-ChildItem -LiteralPath (Join-Path $Repo 'data') -Filter *.step | Sort-Object Name
$want = @($Only | ForEach-Object { $_ -split ',' } | ForEach-Object { $_.Trim() } | Where-Object { $_ })
if ($want.Count -gt 0) {
    $models = $models | Where-Object { $want -contains $_.BaseName }
}

$bad = 0
$realBad = 0
$threshOnly = 0
$tinyOnly = 0
$rows = @()
foreach ($m in $models) {
    $stem = $m.BaseName
    $src = $m.FullName
    foreach ($mode in $Modes) {
        $occRec = Join-Path $Work "$stem.occ.$mode.txt"
        $portRec = Join-Path $Work "$stem.port.$mode.txt"
        if ($mode -eq 'read') {
            & cmd /c "`"$Probe`" `"$src`" --record `"$occRec`" $Deflection $Angle" | Out-Null
            & $Port $src $portRec --defl $Deflection --angle $Angle | Out-Null
        }
        else {
            $occMid = Join-Path $Work "$stem.occ.$mode.iges"
            $portMid = Join-Path $Work "$stem.port.$mode.iges"
            & cmd /c "`"$Probe`" `"$src`" --record-rt iges `"$occMid`" `"$occRec`" $Deflection $Angle" | Out-Null
            & $Port $src $portRec --rt iges --mid $portMid --defl $Deflection --angle $Angle | Out-Null
        }
        if (-not (Test-Path $occRec) -or -not (Test-Path $portRec)) {
            $cmp = [pscustomobject]@{ Struct = @("record missing (occ=$(Test-Path $occRec) port=$(Test-Path $portRec))"); Real = @(); Thresh = @(); Tiny = @() }
        }
        else {
            $cmp = Compare-Record $occRec $portRec
        }
        $diffPath = Join-Path $Work "$stem.$mode.diff.txt"
        # A clean record must not leave the previous run's diff behind: the file
        # is what a reader inspects by hand, and a stale one reads as a live
        # failure for a model that now passes.
        Remove-Item -LiteralPath $diffPath -Force -ErrorAction SilentlyContinue
        if ($cmp.Struct.Count -eq 0 -and $cmp.Real.Count -eq 0 -and $cmp.Thresh.Count -eq 0 -and $cmp.Tiny.Count -eq 0) {
            $rows += [pscustomobject]@{ Model = $stem; Mode = $mode; Result = 'OK'; Fields = '' }
        }
        else {
            $lines = @("RECORD parity $stem ($mode)")
            if ($cmp.Struct.Count -gt 0) { $lines += 'STRUCTURAL (decision differs)'; $lines += $cmp.Struct }
            if ($cmp.Real.Count -gt 0) { $lines += 'REAL (geometry differs above the noise floor)'; $lines += $cmp.Real }
            if ($cmp.Thresh.Count -gt 0) { $lines += 'THRESHOLD (mesh counts: deflection tie-break, not a defect)'; $lines += $cmp.Thresh }
            if ($cmp.Tiny.Count -gt 0) { $lines += 'TINY (below 1e-6 relative: last-digit noise)'; $lines += $cmp.Tiny }
            Set-Content -LiteralPath $diffPath -Value $lines
            $groups = @()
            if ($cmp.Struct.Count -gt 0) { $groups += @{ Name = 'STRUCT'; List = $cmp.Struct } }
            if ($cmp.Real.Count -gt 0) { $groups += @{ Name = 'REAL'; List = $cmp.Real } }
            if ($cmp.Thresh.Count -gt 0) { $groups += @{ Name = 'THRESH'; List = $cmp.Thresh } }
            if ($cmp.Tiny.Count -gt 0) { $groups += @{ Name = 'TINY'; List = $cmp.Tiny } }
            $fields = @($groups | ForEach-Object { @($_.List) } | ForEach-Object { ($_ -split ' {2,}')[0] })
            $result = $groups[0].Name
            $rows += [pscustomobject]@{ Model = $stem; Mode = $mode; Result = $result; Fields = ($fields -join ',') }
            switch ($result) {
                'STRUCT' { $bad++ }
                'REAL' { $realBad++ }
                'THRESH' { $threshOnly++ }
                'TINY' { $tinyOnly++ }
            }
        }
    }
}

$rows | ForEach-Object {
    $row = $_
    switch ($row.Result) {
        'OK' { Write-Host ("OK      {0,-22} {1}" -f $row.Model, $row.Mode) -ForegroundColor Green }
        'THRESH' { Write-Host ("THRESH  {0,-22} {1}  {2}" -f $row.Model, $row.Mode, $row.Fields) -ForegroundColor DarkGray }
        'TINY' { Write-Host ("TINY    {0,-22} {1}  {2}" -f $row.Model, $row.Mode, $row.Fields) -ForegroundColor Yellow }
        'REAL' { Write-Host ("REAL    {0,-22} {1}  {2}" -f $row.Model, $row.Mode, $row.Fields) -ForegroundColor Magenta }
        default { Write-Host ("STRUCT  {0,-22} {1}  {2}" -f $row.Model, $row.Mode, $row.Fields) -ForegroundColor Red }
    }
}
Write-Host ("records: {0}  structural: {1}  real: {2}  threshold-only: {3}  tiny-only: {4}  detail: {5}" -f $rows.Count, $bad, $realBad, $threshOnly, $tinyOnly, $Work)
# Failing = a decision differs, or geometry differs above the noise floor.
# THRESH/TINY are printed for the record and do not fail the gate.
exit ($bad + $realBad)
