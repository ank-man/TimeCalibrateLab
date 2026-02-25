/**
 * Utility functions for sampling and data export.
 */

/**
 * Convert MCMC samples to CSV format for export.
 */
export function samplesToCSV(
  nodeAges: Record<string, number[]>,
  globalRate: number[],
  rates?: Record<string, number[]>
): string {
  const nodeIds = Object.keys(nodeAges);
  const rateIds = rates ? Object.keys(rates) : [];
  const n = globalRate.length;

  const headers = ['iteration', 'globalRate', ...nodeIds.map(id => `age_${id}`), ...rateIds.map(id => `rate_${id}`)];
  const rows = [headers.join(',')];

  for (let i = 0; i < n; i++) {
    const row = [
      i.toString(),
      globalRate[i].toFixed(6),
      ...nodeIds.map(id => (nodeAges[id]?.[i] ?? 0).toFixed(4)),
      ...rateIds.map(id => (rates?.[id]?.[i] ?? 0).toFixed(6)),
    ];
    rows.push(row.join(','));
  }
  return rows.join('\n');
}

/**
 * Export tree as JSON for download.
 */
export function treeToJSON(nodes: Record<string, unknown>): string {
  return JSON.stringify(nodes, null, 2);
}

/**
 * Trigger a file download in the browser.
 */
export function downloadFile(content: string, filename: string, mimeType: string = 'text/plain'): void {
  const blob = new Blob([content], { type: mimeType });
  const url = URL.createObjectURL(blob);
  const link = document.createElement('a');
  link.href = url;
  link.download = filename;
  document.body.appendChild(link);
  link.click();
  document.body.removeChild(link);
  URL.revokeObjectURL(url);
}

/**
 * Download an SVG element as PNG.
 */
export function downloadSVGasPNG(svgElement: SVGSVGElement, filename: string): void {
  const serializer = new XMLSerializer();
  const svgString = serializer.serializeToString(svgElement);
  const canvas = document.createElement('canvas');
  const bbox = svgElement.getBoundingClientRect();
  canvas.width = bbox.width * 2;
  canvas.height = bbox.height * 2;
  const ctx = canvas.getContext('2d');
  if (!ctx) return;

  ctx.scale(2, 2);
  const img = new Image();
  const svgBlob = new Blob([svgString], { type: 'image/svg+xml;charset=utf-8' });
  const url = URL.createObjectURL(svgBlob);

  img.onload = () => {
    ctx.fillStyle = 'white';
    ctx.fillRect(0, 0, canvas.width, canvas.height);
    ctx.drawImage(img, 0, 0);
    URL.revokeObjectURL(url);
    const pngUrl = canvas.toDataURL('image/png');
    const link = document.createElement('a');
    link.href = pngUrl;
    link.download = filename;
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
  };
  img.src = url;
}

/**
 * Format a number for display.
 */
export function formatNum(n: number, decimals: number = 2): string {
  return n.toFixed(decimals);
}

/**
 * Clamp a number between min and max.
 */
export function clamp(value: number, min: number, max: number): number {
  return Math.max(min, Math.min(max, value));
}
