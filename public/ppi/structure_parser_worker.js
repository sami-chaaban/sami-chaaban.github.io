import { parseStructureRecords } from './structure_parser.js';
self.onmessage = event => {
  try { self.postMessage({ data:parseStructureRecords(event.data.text,event.data.format) }); }
  catch (error) { self.postMessage({ error:error?.message || String(error) }); }
};
