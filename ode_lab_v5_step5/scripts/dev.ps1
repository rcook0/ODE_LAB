param(
  [Parameter(Position=0)]
  [ValidateSet("wasm","web","all")]
  [string]$Cmd="all"
)

$Root = Split-Path -Parent $PSScriptRoot
$Web  = Join-Path $Root "web"
$Rust = Join-Path $Root "rust-wasm"

function Wasm {
  Push-Location $Rust
  wasm-pack build --target web --release
  Pop-Location
  if(Test-Path (Join-Path $Web "pkg")) { Remove-Item (Join-Path $Web "pkg") -Recurse -Force }
  Copy-Item (Join-Path $Rust "pkg") (Join-Path $Web "pkg") -Recurse -Force
}

function Web {
  python -m http.server 8080 -d $Web
}

switch($Cmd){
  "wasm" { Wasm }
  "web"  { Web }
  "all"  { Wasm; Web }
}
