// @ts-check

const config = {
  title: 'FedGWAS',
  tagline: 'Federated GWAS screening pipeline documentation',
  url: 'https://idsla.github.io',
  baseUrl: '/Fed-GWAS/',
  organizationName: 'idsla',
  projectName: 'Fed-GWAS',
  onBrokenLinks: 'throw',
  onBrokenMarkdownLinks: 'warn',
  i18n: {
    defaultLocale: 'en',
    locales: ['en'],
  },
  presets: [
    [
      'classic',
      {
        docs: {
          path: 'content/get-started',
          sidebarPath: require.resolve('./sidebarsGetStarted.js'),
          routeBasePath: 'get-started',
          includeCurrentVersion: true,
          lastVersion: 'current',
          versions: {
            current: {
              label: 'v0.0.1',
              path: '',
            },
          },
        },
        blog: false,
        theme: {
          customCss: require.resolve('./src/css/custom.css'),
        },
      },
    ],
  ],
  markdown: {
    mermaid: true,
  },
  plugins: [
    [
      '@docusaurus/plugin-content-docs',
      {
        id: 'userGuide',
        path: 'content/user-guide',
        routeBasePath: 'user-guide',
        sidebarPath: require.resolve('./sidebarsUserGuide.js'),
      },
    ],
    [
      '@docusaurus/plugin-content-docs',
      {
        id: 'examples',
        path: 'content/examples',
        routeBasePath: 'examples',
        sidebarPath: require.resolve('./sidebarsExamples.js'),
      },
    ],
    [
      '@docusaurus/plugin-content-docs',
      {
        id: 'apiReference',
        path: 'content/api-reference',
        routeBasePath: 'api-reference',
        sidebarPath: require.resolve('./sidebarsReference.js'),
      },
    ],
    [
      require.resolve('@easyops-cn/docusaurus-search-local'),
      {
        hashed: true,
        indexDocs: true,
        indexBlog: false,
        indexPages: true,
        language: ['en'],
        docsDir: ['content/get-started', 'content/user-guide', 'content/examples', 'content/api-reference'],
        docsRouteBasePath: ['/get-started', '/user-guide', '/examples', '/api-reference'],
        docsPluginIdForPreferredVersion: 'default',
      },
    ],
  ],
  themes: ['@docusaurus/theme-mermaid'],
  themeConfig: {
    navbar: {
      title: 'FedGWAS',
      items: [
        {
          type: 'docSidebar',
          sidebarId: 'getStartedSidebar',
          docsPluginId: 'default',
          position: 'left',
          label: 'Get Started',
        },
        {
          type: 'docSidebar',
          sidebarId: 'userGuideSidebar',
          docsPluginId: 'userGuide',
          position: 'left',
          label: 'User Guide',
        },
        {
          type: 'docsVersionDropdown',
          docsPluginId: 'default',
          position: 'right',
        },
        {
          type: 'docSidebar',
          sidebarId: 'examplesSidebar',
          docsPluginId: 'examples',
          label: 'Examples',
          position: 'left',
        },
        {
          type: 'docSidebar',
          sidebarId: 'referenceSidebar',
          docsPluginId: 'apiReference',
          label: 'API Reference',
          position: 'left',
        },
      ],
    },
    footer: {
      style: 'dark',
      links: [
        {
          title: 'Project',
          items: [
            {
              label: 'Quick start',
              to: '/get-started/local-simulation',
            },
            {
              label: 'Run the pipeline',
              to: '/user-guide/simulation',
            },
            {
              label: 'Experiments',
              to: '/user-guide/experiments',
            },
            {
              label: 'Examples',
              to: '/examples/overview',
            },
          ],
        },
        {
          title: 'Reference',
          items: [
            {
              label: 'Protocol workflow',
              to: '/user-guide/design/workflow',
            },
            {
              label: 'API reference',
              to: '/api-reference/api',
            },
            {
              label: 'Output files',
              to: '/api-reference/outputs',
            },
            {
              label: 'Flower',
              href: 'https://flower.ai',
            },
            {
              label: 'PLINK',
              href: 'https://www.cog-genomics.org/plink/1.9/',
            },
          ],
        },
      ],
      copyright: `Copyright (c) ${new Date().getFullYear()} FedGWAS contributors.`,
    },
    prism: {
      additionalLanguages: ['bash', 'python', 'yaml', 'toml'],
    },
    colorMode: {
      defaultMode: 'light',
      respectPrefersColorScheme: true,
    },
  },
};

module.exports = config;
