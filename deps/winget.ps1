$url   = $args[0]
$fname = $args[1]
[Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12 

(new-object net.webclient).DownloadFile($url, $fname)
