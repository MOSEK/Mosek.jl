$url   = $args[0]
$fname = $args[1]

(new-object net.webclient).DownloadFile($url, $filename)
