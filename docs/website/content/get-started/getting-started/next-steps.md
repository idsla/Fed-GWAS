---
title: Next Steps
slug: /next-steps
---

# Next Steps

After you have read the overview and completed either a local simulation or a
federated deployment setup, use the rest of the documentation based on what you
want to do next.

## Configure a study

Use the User Guide when you are preparing your own data, adjusting runtime
settings, or operating FedGWAS beyond the first quick run.

- [Configuration](/user-guide/configuration): center YAML files, input data paths,
  thresholds, output directories, and Flower settings.
- [Simulation Mode](/user-guide/simulation): repository-level local simulation behavior,
  Flower run config, scenario layout, and output locations.
- [Installation and Overview](/user-guide/federated-deployment/installation):
  federated roles, ports, and checks.
- [CLI Deployment](/user-guide/federated-deployment/cli-deployment) and
  [Script Deployment](/user-guide/federated-deployment/script-deployment):
  the two ways to start SuperLink, SuperNodes, and `flwr run`.
- [fedgwas-sim CLI](/user-guide/cli-simulation): full command reference for local
  simulation projects, presets, checks, evaluation, and result collection.
- [Troubleshooting](/user-guide/troubleshooting): common setup, data, PLINK, and runtime
  issues.

## Run examples

Use Examples when you want a known scenario before moving to your own study.

- [Examples overview](/examples/overview): choose between correctness,
  performance, and real-world scenarios.
- [Tiny correctness](/examples/tiny-correctness): minimal synthetic run for
  validating the end-to-end pipeline.
- [Performance small](/examples/performance-small): small benchmark-style run
  for timing and resource checks.
- [1000 Genomes](/examples/1000genomes): real-world chromosome 22 workflow.

## Understand the design

Use Design when you need to understand why the pipeline is structured the way it
is, or when you are reviewing the protocol before extending the code.

- [Architecture](/user-guide/design/architecture): major components and how Flower,
  clients, server strategy, PLINK, and output directories fit together.
- [Workflow](/user-guide/design/workflow): four conceptual screening stages and
  the federated rounds from initialization through association screening.
- [Privacy masking](/user-guide/design/privacy-masking): encryption, shuffling,
  anonymization, lightweight secret-sharing, and the server’s relay role.

## Reference APIs and outputs

Use API Reference when you need exact runtime fields, module responsibilities,
or output file expectations.

- [API Reference](/api-reference/api): Flower run config, server and center
  configuration fields, and main Python entry points.
- [Output files](/api-reference/outputs): result directories and generated
  artifacts.
- [Module Reference](/api-reference/modules/client): client, server, quality
  control, kinship, and association modules.
