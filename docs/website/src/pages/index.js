import React from 'react';
import clsx from 'clsx';
import Link from '@docusaurus/Link';
import useDocusaurusContext from '@docusaurus/useDocusaurusContext';
import Layout from '@theme/Layout';
import Heading from '@theme/Heading';
import styles from './index.module.css';

const features = [
  {
    title: 'Federated GWAS screening',
    description:
      'Coordinate federated quality control, KING-based relatedness screening, and association screening across participating clients so that each site can shrink the set of loci that still need follow-up testing.',
  },
  {
    title: 'Privacy-preserving coordination',
    description:
      'Keep genotype-level computation local while using encryption, shuffling, anonymization, and lightweight secret-sharing. The server relays selected protocol messages without decrypting them.',
  },
  {
    title: 'Experiment-ready layout',
    description:
      'Run repeatable Flower simulation and local-deployment experiments with scenario-specific configs, logs, and result directories.',
  },
];

const capabilities = [
  {
    title: 'Encrypted relay',
    description:
      'Clients exchange encrypted seed shares and other protocol payloads via the server, which forwards messages without decrypting them.',
    accent: 'Server as relay',
  },
  {
    title: 'Federated QC',
    description:
      'Harmonize sample missingness, SNP missingness, minor allele frequency, and Hardy–Weinberg equilibrium filters across participating clients.',
    accent: 'PLINK-backed QC',
  },
  {
    title: 'Relatedness Screening',
    description:
      'Estimate pairwise kinship with KING and optionally remove related samples before association screening.',
    accent: 'KING kinship',
  },
  {
    title: 'Association Screening',
    description:
      'Run local logistic-regression filtering with privacy-preserving tokens, then federated case–control logistic regression, while keeping raw genotypes at each client.',
    accent: 'Locus prioritization',
  },
];

const highlights = ['Local data ownership', 'Privacy-preserving relay', 'GWAS screening'];

const entryPoints = [
  {
    title: 'Prerequisites',
    description: 'Set up Python, PLINK, Flower, and the required runtime tools.',
    to: '/get-started/prerequisites',
  },
  {
    title: 'Local Simulation',
    description: 'Run the default tiny correctness experiment from data generation to evaluation.',
    to: '/get-started/local-simulation',
  },
  {
    title: 'Examples',
    description: 'Browse tiny correctness, performance small, and 1000 Genomes workflows.',
    to: '/examples/overview',
  },
  {
    title: 'API Reference',
    description: 'Look up Flower run config, YAML schema, output files, and Python entry points.',
    to: '/api-reference/api',
  },
];

function HomepageHeader() {
  const {siteConfig} = useDocusaurusContext();
  return (
    <header className={clsx('hero', styles.hero)}>
      <div className="container">
        <div className={styles.heroGrid}>
          <div className={styles.heroCopy}>
            <Heading as="h1" className={styles.heroTitle}>
              {siteConfig.title}
            </Heading>
            <p className={styles.heroSubtitle}>
              A lightweight federated pipeline for
              privacy-preserving GWAS screening across distributed genomic
              datasets. Clients keep genotype-level data local while coordinating
              federated quality control, relatedness screening, and association screening.
            </p>
            <div className={styles.highlights} aria-label="Project highlights">
              {highlights.map((highlight) => (
                <span key={highlight}>{highlight}</span>
              ))}
            </div>
            <div className={styles.actions}>
              <Link className="button button--primary button--lg" to="/get-started/">
                Get Started
              </Link>
              <Link className="button button--secondary button--lg" to="/api-reference/api">
                API reference
              </Link>
            </div>
          </div>
          <div className={styles.pipelinePanel} aria-label="Pipeline capabilities">
            <div className={styles.panelHeader}>
              <span>Screening Capabilities</span>
              <span className={styles.status}>Privacy-preserving</span>
            </div>
            <div className={styles.capabilityGrid}>
              {capabilities.map((capability) => (
                <div className={styles.capability} key={capability.title}>
                  <div>
                    <span className={styles.capabilityAccent}>{capability.accent}</span>
                    <h2>{capability.title}</h2>
                    <p>{capability.description}</p>
                  </div>
                </div>
              ))}
            </div>
          </div>
        </div>
      </div>
    </header>
  );
}

function FeatureGrid() {
  return (
    <section className={styles.features}>
      <div className="container">
        <div className={styles.sectionHeader}>
          <Heading as="h2">Built for research operations</Heading>
          <p>
            The documentation follows the repository structure so researchers can move from setup
            to pipeline internals without switching mental models.
          </p>
        </div>
        <div className={styles.cards}>
          {features.map((feature) => (
            <article className={styles.card} key={feature.title}>
              <h3>{feature.title}</h3>
              <p>{feature.description}</p>
            </article>
          ))}
        </div>
      </div>
    </section>
  );
}

function EntryPoints() {
  return (
    <section className={styles.entryPoints}>
      <div className="container">
        <div className={styles.sectionHeader}>
          <Heading as="h2">Start from the task in front of you</Heading>
          <p>
            Core documentation paths are split by workflow so setup, examples, and reference
            material stay easy to scan.
          </p>
        </div>
        <div className={styles.entryGrid}>
          {entryPoints.map((entry) => (
            <Link className={styles.entryCard} to={entry.to} key={entry.title}>
              <h3>{entry.title}</h3>
              <p>{entry.description}</p>
            </Link>
          ))}
        </div>
      </div>
    </section>
  );
}

function QuickLinks() {
  return (
    <section className={styles.quickLinks}>
      <div className="container">
        <div className={styles.quickGrid}>
          <Link to="/user-guide/configuration">Configure centers</Link>
          <Link to="/user-guide/design/workflow">Review workflow</Link>
          <Link to="/api-reference/outputs">Interpret outputs</Link>
          <Link to="/user-guide/troubleshooting">Troubleshoot runs</Link>
        </div>
      </div>
    </section>
  );
}

export default function Home() {
  const {siteConfig} = useDocusaurusContext();
  return (
    <Layout
      title={`${siteConfig.title} documentation`}
      description="A lightweight federated pipeline for privacy-preserving GWAS screening across distributed genomic datasets">
      <HomepageHeader />
      <main>
        <EntryPoints />
        <FeatureGrid />
        <QuickLinks />
      </main>
    </Layout>
  );
}
