// @ts-check
import { defineConfig } from 'astro/config';

function normalizeBase(value) {
  if (!value || value === '/') return '/';
  let base = value;
  if (!base.startsWith('/')) base = `/${base}`;
  if (base.endsWith('/')) base = base.slice(0, -1);
  return base;
}

const repo = process.env.GITHUB_REPOSITORY?.split('/')[1] ?? '';
const owner = process.env.GITHUB_REPOSITORY?.split('/')[0] ?? '';
const isUserSite = owner && repo && repo === `${owner}.github.io`;

const defaultBase = isUserSite ? '/' : repo ? `/${repo}` : '/';
const defaultSite =
  owner && repo
    ? `https://${owner}.github.io${isUserSite ? '' : `/${repo}`}`
    : '';

const site = process.env.SITE_URL || defaultSite || undefined;
const base = normalizeBase(process.env.BASE_PATH || defaultBase);

// https://astro.build/config
export default defineConfig({
  site,
  base,
  output: 'static',
});
