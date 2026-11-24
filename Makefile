IMAGE ?= 3dtrees_trajectory:dev
TEST_IMAGE ?= $(IMAGE)-test
TEST_CPUS ?= 2
TEST_MEMORY ?= 4g
TEST_CREATED_BY ?= trajectory-tests
DOCKER_BUILD_FLAGS ?=
DOCKER_RUN_FLAGS ?=

.PHONY: image test
image:
	docker build $(DOCKER_BUILD_FLAGS) --target runtime --tag $(IMAGE) .

test:
	docker build $(DOCKER_BUILD_FLAGS) --target test --tag $(TEST_IMAGE) .
	docker run --rm --name trajectory-tests-$$$$ --label "owner=$$(id -un)" \
		--label "created-by=$(TEST_CREATED_BY)" --network none --cpus $(TEST_CPUS) \
		--memory $(TEST_MEMORY) --memory-swap $(TEST_MEMORY) \
		-e OPENBLAS_NUM_THREADS=1 -e OMP_NUM_THREADS=1 \
		-e MKL_NUM_THREADS=1 -e RAYON_NUM_THREADS=$(TEST_CPUS) \
		$(DOCKER_RUN_FLAGS) $(TEST_IMAGE)
