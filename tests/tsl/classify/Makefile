ROOT = ../../..

# force sort to be consistent
export LC_ALL=C


binDir = ${ROOT}/bin

test:: unitTests

unitTests:
	python classifyUnitTests.py


mkdirs:
	@mkdir -p output

clean::
	rm -rf output
