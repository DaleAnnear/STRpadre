# Test strategy

Run locally with pinned test dependencies:

```bash
python -m pip install -r requirements-test.txt
pytest -q
ruff check bin tests
python -m py_compile bin/*.py
```

The unit suite covers schema/configuration rejection, bad caller selection, native VCF
parsers, haploid and filtered calls, min-cost unordered diploid pairing, motif rotations
and reverse complements, high discordance, and mixed-platform caller eligibility.

The Nextflow stub test exercises caller selection and output wiring without executing a
bioinformatics executable:

```bash
bash tests/test_nextflow_stub.sh
```

The small Docker integration builds only the normalizer/consensus image and executes a
real normalizer plus consensus test over public synthetic VCF fixtures:

```bash
bash tests/test_docker_integration.sh
```

Caller executable smoke tests require images built under their license conditions:

```bash
bash containers/smoke_test.sh
```
