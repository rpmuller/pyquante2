.PHONY: build test
build:
	uv pip install -e ".[dev]"

test:
	uv run pytest --doctest-modules --ignore=othertests
