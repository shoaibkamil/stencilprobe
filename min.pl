$a = 999999.9999;
while (<>) {
	if ($_ < $a) {
		$a = $_;
	}
}
print $a;
