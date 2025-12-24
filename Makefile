build:
	pip install -e ".[dev]"

test:
	pytest --doctest-modules --ignore=othertests
