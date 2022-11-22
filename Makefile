test:
	pytest

test-cov:
	pytest --cov=cov_pairs --cov-report=html

clean:
	rm -rf htmlcov && rm .coverage
