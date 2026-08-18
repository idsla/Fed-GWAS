module.exports = {
  getStartedSidebar: [
    {
      type: 'doc',
      id: 'intro',
      label: 'FedGWAS overview',
    },
    {
      type: 'doc',
      id: 'getting-started/prerequisites',
      label: 'Prerequisites',
    },
    {
      type: 'doc',
      id: 'getting-started/local-simulation',
      label: 'Local Simulation',
    },
    {
      type: 'doc',
      id: 'getting-started/federated-deployment',
      label: 'Federated Deployment',
    },
    {
      type: 'doc',
      id: 'getting-started/next-steps',
      label: 'Next Steps',
    },
  ],
  userGuideSidebar: [
    {
      type: 'category',
      label: 'User Guide',
      collapsed: false,
      items: [
        'getting-started/configuration',
        'getting-started/simulation',
        'cli_simulation',
        'getting-started/deployment',
      ],
    },
    {
      type: 'category',
      label: 'Experiments and Operations',
      collapsed: false,
      items: [
        'experiments/overview',
        'experiments/runner',
        'troubleshooting',
      ],
    },
    {
      type: 'category',
      label: 'Design',
      collapsed: false,
      items: [
        'design/architecture',
        'design/workflow',
        'design/privacy-masking',
      ],
    },
  ],
};
