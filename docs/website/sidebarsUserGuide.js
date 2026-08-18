module.exports = {
  userGuideSidebar: [
    {
      type: 'category',
      label: 'Local Simulation',
      collapsed: false,
      items: [
        'local-simulation/installation',
        'local-simulation/config-and-data',
        'local-simulation/run-experiment',
        'local-simulation/evaluation',
        'local-simulation/cli_simulation',
      ],
    },
    {
      type: 'category',
      label: 'Federated Deployment',
      collapsed: false,
      items: [
        {
          type: 'doc',
          id: 'federated-deployment/installation',
          label: 'Installation and Overview',
        },
        {
          type: 'doc',
          id: 'federated-deployment/cli-deployment',
          label: 'CLI Deployment',
        },
        {
          type: 'doc',
          id: 'federated-deployment/script-deployment',
          label: 'Script Deployment',
        },
      ],
    },
    {
      type: 'category',
      label: 'Design',
      collapsed: false,
      items: [
        'design/architecture',
        'design/workflow',
        'design/gwas-components',
        'design/parameters',
        'design/privacy-masking',
      ],
    },
    {
      type: 'doc',
      id: 'troubleshooting',
      label: 'Troubleshooting',
    },
  ],
};
