build:
	docker build -t uscbiostats/fmcmc-dev:latest --build-arg Rversion=latest -f Dockerfile.r-ver . ; \
		docker build -t uscbiostats/fmcmc-dev:devel -f Dockerfile.r-dev .
push: build
	docker push uscbiostats/fmcmc-dev:latest ; \
		docker push uscbiostats/fmcmc-dev:devel
