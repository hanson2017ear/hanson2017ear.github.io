/* TCGA LIHC exploratory reader. All model estimates come from the supplied R output.
   The browser only selects, filters, scales and displays those estimates. */
(() => {
  'use strict';
  const D = window.APP_DATA;
  // Only presentation copy is merged. Numerical results and original R excerpts are unchanged.
  function mergeCopy(target, patch) {
    if (!target || !patch) return;
    for (const key of Object.keys(patch)) {
      if (patch[key] && typeof patch[key] === 'object') mergeCopy(target[key], patch[key]);
      else target[key] = patch[key];
    }
  }
  mergeCopy(D, window.SITE_COPY);
  const $ = (id) => document.getElementById(id);
  const esc = (s) => String(s ?? '').replace(/[&<>"']/g, (c) => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
  const finite = (x) => x !== null && x !== undefined && Number.isFinite(Number(x));
  const fmt = (x, digits = 2) => finite(x) ? Number(x).toLocaleString('en-US', {maximumFractionDigits:digits}) : 'NA';
  const sci = (x) => !finite(x) ? 'NA' : Number(x) === 0 ? '0' : Number(x) < .001 ? Number(x).toExponential(2) : Number(x).toPrecision(3);
  const sign = (x) => finite(x) ? `${Number(x) > 0 ? '+' : ''}${Number(x).toFixed(3)}` : 'NA';
  const colors = {GroupX:'#D55E00',GroupY:'#0072B2',GroupX_higher:'#D55E00',GroupY_higher:'#0072B2'};
  const ids = ['overview','quality','design','samples','differential','heatmap','ora','gsea','interpretation'];
  const navNames = ['研究概览','数据检查','配对模型','样本结构','差异表达','基因热图','通路富集','全基因 GSEA','证据与判定'];
  const topics = ['先理解问题与数据','原始数据检查','配对设计','样本结构','差异表达','基因热图','GO / KEGG','全基因 GSEA','组织类型判定'];
  let active = 0, codeMode = 'step', currentCode = '', volcano = null, enrichmentState = null, redraw = null, canvasObserver = null;
  const state = {gene:'EZH2',fdr:.05,lfc:1,patient:null,sampleTab:'pca',heatOrder:'group'};
  const runCode = '# 下载本页的 TCGA_LIHC_analysis.R\n# 输入请使用课件原附的两份 CSV：\n# TCGA_LIHC的8个样本count.csv\n# TCGA_LIHC的8个样本分组.csv\n# 本网站不提供原始输入档下载。\n\n# 在 R 中先依完整代码开头的说明设定工作目录。\n# 代码会检查依赖，并使用原始 count 执行配对分析。\nsource("TCGA_LIHC_analysis.R", encoding = "UTF-8")\n\n# 正式比较方向\n# log2FoldChange > 0：GroupX 较高\n# NES > 0：GroupX 方向富集\n# NES < 0：GroupY 方向富集\n\n# 网站互动只使用既有分析结果。\n# 调整阈值不会重新拟合模型或重新做 BH 校正。';

  if (!D || !Array.isArray(D.steps) || !Array.isArray(D.de)) {
    $('step-content').innerHTML = '<div class="error-state"><h1>分析数据尚未加载</h1><p>请确认 app-data.js 与本页位于同一数据夹，再刷新。网站需要正式分析数据才能显示结果。</p></div>';
    $('previous').disabled = $('next').disabled = true;
    return;
  }
  const counts = D.meta?.counts || {};
  const de = D.de;
  const genes = new Map(de.map((g) => [String(g.gene).toUpperCase(),g]));
  const samples = D.samples || [];
  const getSample = (id) => samples.find((s) => [s.id,s.label,s.sample].includes(id));
  const groupOf = (id) => getSample(id)?.group || (String(id).endsWith('_X') ? 'GroupX' : 'GroupY');
  const patientOf = (id) => getSample(id)?.patient || String(id).split('_')[0];
  const countVal = (key, fallback) => finite(counts[key]) ? Number(counts[key]) : fallback;
  const cx = countVal('degX',de.filter((g) => finite(g.padj) && g.padj < .05 && g.lfc > 1).length);
  const cy = countVal('degY',de.filter((g) => finite(g.padj) && g.padj < .05 && g.lfc < -1).length);
  const figures = D.figures || D.steps.flatMap((s) => s.figures || []);
  const renderFindings = (step) => (step.findings?.length ? `<section class="result-reading"><h2>本步结果</h2><ul class="findings">${step.findings.map((f) => `<li>${esc(f)}</li>`).join('')}</ul></section>` : '');
  const figureButtons = (fs, heading = '查看原始图表') => !fs?.length ? '' : `<section class="figure-section"><h3>${heading}</h3><div class="originals">${fs.map((f) => `<button class="figure-button" data-figure="${esc(f.png)}" data-title="${esc(f.title)}" data-pdf="${esc(f.pdf || '')}">${esc(f.title)}<span aria-hidden="true">↗</span></button>`).join('')}</div></section>`;
  const fdrTag = (fdr) => `<span title="BH 校正后 P 值">${sci(fdr)}</span>`;
  const groupTag = (d) => `<span class="${d === 'GroupX_higher' ? 'x-text' : 'y-text'}">${d === 'GroupX_higher' ? 'GroupX' : 'GroupY'}</span>`;

  function overviewPca() {
    const points=D.pca.points,w=480,h=316,p={l:37,r:23,t:25,b:43};
    const xx=(v)=>p.l+(v+56)/204*(w-p.l-p.r), yy=(v)=>h-p.b-(v+65)/108*(h-p.t-p.b);
    let out=`<svg viewBox="0 0 ${w} ${h}" role="img" aria-label="八个真实样本的 PCA。灰线连接同一患者，橙色为 GroupX，蓝色为 GroupY。">`;
    for (const v of [-50,0,50,100]) out+=`<line x1="${xx(v)}" y1="${p.t}" x2="${xx(v)}" y2="${h-p.b}" stroke="#dde4d4" stroke-dasharray="2 5"/><text x="${xx(v)}" y="${h-p.b+18}" text-anchor="middle" fill="#697a61" font-size="9">${v}</text>`;
    for (const v of [-50,-25,0,25]) out+=`<line x1="${p.l}" y1="${yy(v)}" x2="${w-p.r}" y2="${yy(v)}" stroke="#dde4d4" stroke-dasharray="2 5"/><text x="${p.l-10}" y="${yy(v)+3}" text-anchor="end" fill="#697a61" font-size="9">${v}</text>`;
    for(const patient of ['P1','P2','P3','P4']){const pair=points.filter(v=>v.patient===patient);out+=`<line x1="${xx(pair[0].pc1)}" y1="${yy(pair[0].pc2)}" x2="${xx(pair[1].pc1)}" y2="${yy(pair[1].pc2)}" stroke="#b5c4a9" stroke-width="1.1"/>`;}
    for(const v of points){const x=xx(v.pc1),y=yy(v.pc2);out+=`<circle cx="${x}" cy="${y}" r="5" fill="${colors[v.group]}" stroke="#f4f6ee" stroke-width="2"/>`;if(v.group==='GroupX')out+=`<text x="${x+(v.patient==='P4'?-11:10)}" y="${y+(v.patient==='P3'?19:-10)}" text-anchor="${v.patient==='P4'?'end':'start'}" fill="#9c511b" font-size="10">${esc(v.sample)}</text>`;}
    out+=`<text x="${xx(-38)}" y="${yy(33)}" text-anchor="middle" fill="#266887" font-size="10">GroupY</text><text x="${(w+p.l-p.r)/2}" y="${h-4}" text-anchor="middle" fill="#66795c" font-size="9">PC1 · ${fmt(D.pca.variance[0],2)}%</text><text transform="translate(11 ${(h-p.b+p.t)/2}) rotate(-90)" text-anchor="middle" fill="#66795c" font-size="9">PC2 · ${fmt(D.pca.variance[1],2)}%</text></svg>`;
    return out;
  }
  function overview() {
    const patients=[...new Set(samples.map(s=>s.patient))].sort();
    return `<section class="overview-hero"><div class="hero-copy"><div class="hero-eyebrow">TCGA–LIHC / RNA-seq</div><h1>从八个样本，<br><em>判断肿瘤与癌旁。</em></h1><p class="hero-description">数据来自 4 位患者，每位患者各有 GroupX 和 GroupY 一例样本。先比较表达和通路，再核对组织类型。</p><div class="hero-actions"><a class="button" href="#quality">开始阅读分析 <span aria-hidden="true">→</span></a><a class="hero-link" href="#interpretation">查看组别判定 <span aria-hidden="true">↗</span></a></div><p class="hero-method">DESeq2 · GO / KEGG · GSEA</p></div><div class="hero-plot"><div class="hero-plot-head"><strong>八个样本的表达结构</strong><span>VST / PCA</span></div><a href="#samples" aria-label="查看 PCA 详细结果与患者配对">${overviewPca()}</a><div class="hero-plot-bottom"><div class="group-key"><span class="key-x">GroupX</span><span class="key-y">GroupY</span></div><a href="#samples">查看配对样本 ↗</a></div></div></section>
    <div class="facts"><div class="fact"><span>配对样本</span><strong>${fmt(countVal('samples',samples.length),0)}<small>个</small></strong></div><div class="fact"><span>患者配对</span><strong>${fmt(patients.length,0)}<small>对</small></strong></div><div class="fact"><span>进入模型的基因</span><strong>${fmt(de.length,0)}</strong></div><div class="fact"><span>显著差异基因</span><strong>${fmt(cx+cy,0)}</strong></div></div>
    <section class="reading-path"><h2>分析顺序</h2><div><p>先检查原始 count，并确认患者配对；再用 DESeq2 比较两组表达，通过 GO、KEGG 和 GSEA 检查通路。每一步都可以查看图表、筛选结果和对应的 R 代码。</p><a href="#design">查看配对模型 ↗</a></div></section>
    <section class="section"><div class="section-head"><h2><span class="section-number">01</span>同一患者，两种组织</h2><small>每位患者各提供一例 X 和 Y 样本</small></div><div class="paired-grid">${patients.map(p=>`<div class="pair-column"><div class="pair-name">患者 ${esc(p)}</div><div class="pair-sample"><i></i>${esc(p)}_X</div><div class="pair-connector"></div><div class="pair-sample y"><i></i>${esc(p)}_Y</div></div>`).join('')}</div></section><p class="provenance">课程提供的 TCGA–LIHC 8 样本数据 · 正式分析日期 ${esc(D.meta.analysisDate)}<br>页面以 P1–P4 标记患者。统计计算在 R 中完成，网页筛选只调整显示。</p>`;
  }
  function quality() {
    const q = D.qc || {};
    const raw = countVal('inputGenes',q.input_genes);
    const fractional = countVal('excludedFractional',q.fractional_genes);
    const low = countVal('lowExpression',q.low_expression_excluded);
    const modeled = countVal('modeledGenes',de.length);
    const rows = [['样本名称与顺序','矩阵栏名与分组表一一对应，读入后再次依矩阵排序。'],['负数与缺失值',`${fmt(q.negative ?? 0,0)} 个负值；${fmt(q.missing ?? 0,0)} 个缺失值。`],['基因 ID','基因符号形式（SYMBOL），无重复或空白 ID。受检基因的映射结果见富集步骤。'],['真正小数',`${fmt(q.fractional_cells ?? 470,0)} 个值涉及 ${fmt(fractional,0)} 个基因；主分析排除这些基因。`],['微小浮点尾数',`${fmt(q.near_integer_cells ?? 630,0)} 个值在 10⁻⁶ 容差内修为整数。`],['低表达过滤','至少 4 / 8 个样本的 count ≥ 10。'],['原始档案','保留原始文件与杂凑；所有处理均另存结果。']];
    return `<section class="section"><h2>从 ${fmt(raw,0)} 到 ${fmt(modeled,0)} 个基因</h2><div class="retention-track" aria-label="基因过滤比例"><i style="width:${fractional/raw*100}%;background:#bf937b"></i><i style="width:${low/raw*100}%;background:#d9e1e7"></i><i style="width:${modeled/raw*100}%;background:#3c617e"></i></div><div class="retention-labels"><div><i style="background:#bf937b"></i>排除含真正小数的基因<strong>${fmt(fractional,0)}</strong></div><div><i style="background:#d9e1e7"></i>排除低表达基因<strong>${fmt(low,0)}</strong></div><div><i style="background:#3c617e"></i>保留有效整数基因建模<strong>${fmt(modeled,0)}</strong></div></div></section><section class="section"><h2>输入检查纪录</h2><div class="table-scroll"><table class="check-table"><thead><tr><th>检查项目</th><th>正式分析处理</th></tr></thead><tbody>${rows.map((r) => `<tr><td>${r[0]}</td><td>${r[1]}</td></tr>`).join('')}</tbody></table></div></section>`;
  }
  function design() {
    const patients = [...new Set(samples.map((s) => s.patient))].sort();
    return `<div class="design-formula">~ patient + Group</div><p class="plot-hint">以患者作为配对区组，在同一个负二项模型中估计 Group 效应。</p><section class="section"><h2>配对关系</h2><div class="table-scroll"><table class="paired-table"><thead><tr><th>患者</th><th class="x-text">GroupX</th><th class="y-text">GroupY</th></tr></thead><tbody>${patients.map((p) => `<tr><td>${esc(p)}</td><td>${esc(samples.find((s) => s.patient === p && s.group === 'GroupX')?.label || `${p}_X`)}</td><td>${esc(samples.find((s) => s.patient === p && s.group === 'GroupY')?.label || `${p}_Y`)}</td></tr>`).join('')}</tbody></table></div></section><section class="section"><h2>固定比较方向</h2><p class="evidence-text"><span class="inline-code">contrast = c("Group", "GroupX", "GroupY")</span></p><div class="evidence-split"><div><div class="direction-label x-text">log2FoldChange &gt; 0</div><h2>GroupX 表达较高</h2><p>正的 Wald statistic 位于 GSEA 排名左侧，正 NES 对应 GroupX 方向。</p></div><div><div class="direction-label y-text">log2FoldChange &lt; 0</div><h2>GroupY 表达较高</h2><p>负的 Wald statistic 排在排名右侧，负 NES 对应 GroupY 方向。</p></div></div><p class="evidence-text">主模型使用原始整数 count。VST 只用于样本结构和热图，不能代替差异模型的输入。</p></section>`;
  }
  function sampleStructure() {
    return `<div class="controls"><div class="segmented" aria-label="样本图类型"><button data-sample-tab="pca" aria-pressed="${state.sampleTab === 'pca'}">PCA</button><button data-sample-tab="correlation" aria-pressed="${state.sampleTab === 'correlation'}">样本相关性</button></div><div class="patient-buttons" aria-label="选择患者配对"><button data-patient="" aria-pressed="${!state.patient}">全部</button>${['P1','P2','P3','P4'].map((p) => `<button data-patient="${p}" aria-pressed="${state.patient === p}">${p}</button>`).join('')}</div></div><div class="plot-shell"><div class="plot-heading"><h3 id="sample-plot-title"></h3><small>VST 转换后数据</small></div><div id="sample-plot" class="plot-svg-scroll"></div><div id="sample-readout" class="sample-readout">点击患者或图上的样本，追踪同一位患者的配对。</div></div><p class="plot-hint">PCA 的位置来自正式 R 输出。选取患者只改变强调显示，不会重新计算座标。</p>`;
  }
  function differential() {
    return `<div class="controls"><div class="search-wrap"><input type="search" id="gene-search" placeholder="搜索基因，例如 EZH2、FOXM1、GPC3" aria-label="搜索基因" autocomplete="off"><div class="suggestions" id="gene-suggestions" hidden></div></div><button class="figure-button" id="reset-thresholds">回到正式阈值</button></div><div class="plot-shell"><div class="plot-heading"><h3>GroupX 相对 GroupY 的火山图</h3><small>${fmt(de.length,0)} 个受检基因</small></div><div class="sliders"><label><span class="slider-top"><span>padj 小于</span><output id="fdr-value">0.050</output></span><input id="fdr-slider" aria-label="探索 padj 阈值" type="range" min="0.001" max="0.100" step="0.001" value="${state.fdr}"></label><label><span class="slider-top"><span>|log2FC| 大于</span><output id="lfc-value">1.0</output></span><input id="lfc-slider" aria-label="探索 log2FC 阈值" type="range" min="0" max="4" step="0.1" value="${state.lfc}"></label></div><canvas id="volcano" aria-label="可探索的火山图，使用上方搜索栏也可以查询任一受检基因"></canvas><div class="plot-legend"><span class="key-x">GroupX 较高 <strong id="count-x"></strong></span><span class="key-y">GroupY 较高 <strong id="count-y"></strong></span><span class="neutral">其余 <strong id="count-other"></strong></span></div></div><p class="exploration-note"><strong id="threshold-state">正式阈值</strong>　正式结果固定为 padj &lt; 0.05 且 |log2FC| &gt; 1，GroupX ${fmt(cx,0)}、GroupY ${fmt(cy,0)} 个基因。滑杆仅调整显示与计数，不重做模型或 BH 校正。</p><div class="gene-inspector" id="gene-inspector"></div><p class="plot-hint">移到点上看数值；点击可固定基因。也可以在搜索栏输入基因名称后按 Enter。图上不显示 NA 校正值。</p><section class="section"><h2>从正式差异表查阅基因</h2><div class="table-scroll"><table><thead><tr><th>基因</th><th class="numeric">log2FC</th><th class="numeric">BH padj</th><th>较高组别</th></tr></thead><tbody>${de.filter((g) => finite(g.padj) && g.padj < .05 && Math.abs(g.lfc) > 1).sort((a,b) => a.padj-b.padj).slice(0,8).map((g) => `<tr><td><button class="table-row-button" data-gene="${esc(g.gene)}">${esc(g.gene)}</button></td><td class="numeric">${sign(g.lfc)}</td><td class="numeric">${sci(g.padj)}</td><td>${groupTag(g.lfc>0?'GroupX_higher':'GroupY_higher')}</td></tr>`).join('')}</tbody></table></div></section>`;
  }
  function heatmap() {
    return `<div class="controls"><div class="segmented" aria-label="热图样本排序"><button data-heat-order="group" aria-pressed="${state.heatOrder === 'group'}">按组别排列</button><button data-heat-order="patient" aria-pressed="${state.heatOrder === 'patient'}">按患者配对</button></div><span class="plot-hint">移到格子上查看基因与 Z-score</span></div><div class="plot-shell"><div class="plot-heading"><h3>前 40 个显著基因</h3><small>VST，逐基因 Z-score</small></div><div class="plot-svg-scroll" id="heatmap-plot"></div></div><div class="heat-legend"><i></i>相对低表达<i class="z0"></i>0<i class="zplus"></i>相对高表达</div><p class="exploration-note">互动版依数据表顺序展示 40 个基因，仅切换样本列的排列。正式群聚热图可由下方开启。每一列独立标准化，色彩不能用来比较不同基因的绝对表达量。</p>`;
  }
  function enrichment(kind) {
    enrichmentState = enrichmentState?.kind === kind ? enrichmentState : {kind,db:'GO_BP',direction:'GroupX_higher',sig:true,query:'',page:0};
    return `<div class="controls"><div class="segmented" aria-label="富集数据库"><button data-db="GO_BP" aria-pressed="${enrichmentState.db === 'GO_BP'}">GO BP</button><button data-db="KEGG" aria-pressed="${enrichmentState.db === 'KEGG'}">KEGG</button></div><div class="segmented" aria-label="富集方向"><button class="x" data-direction="GroupX_higher" aria-pressed="${enrichmentState.direction === 'GroupX_higher'}">GroupX</button><button class="y" data-direction="GroupY_higher" aria-pressed="${enrichmentState.direction === 'GroupY_higher'}">GroupY</button><button data-direction="all" aria-pressed="${enrichmentState.direction === 'all'}">两组</button></div></div><div id="enrichment-chart"></div><div class="controls"><div class="search-wrap"><input id="term-search" type="search" aria-label="搜索通路名称或 ID" placeholder="搜索通路名称或 ID" value="${esc(enrichmentState.query)}"></div><label><input type="checkbox" id="sig-only" ${enrichmentState.sig?'checked':''}>只看 BH &lt; 0.05</label></div><div class="status-line"><span id="enrichment-status"></span><span>${kind === 'gsea'?'NES 正：GroupX；负：GroupY':'分方向的超几何检验'}</span></div><div class="table-scroll"><table id="enrichment-table"></table></div><div class="table-footer"><span>校正后 P 值来自正式检验</span><div class="pager"><button id="table-prev">上一页</button><span id="table-page"></span><button id="table-next">下一页</button></div></div><div id="gsea-curves"></div>`;
  }
  function interpretation() {
    const wanted = ['EZH2','FOXM1','GPC3','BIRC5','MKI67','CYP2C19','CYP3A4','BCHE'].map((g) => genes.get(g)).filter(Boolean);
    return `<section class="evidence-split"><div><div class="direction-label x-text">GroupX</div><h2>肿瘤组织</h2><p>细胞周期与染色体分离相关表达较高，多个促癌或增殖相关基因共同支持判定。</p></div><div><div class="direction-label y-text">GroupY</div><h2>癌旁非肿瘤组织</h2><p>肝脏代谢程序、补体及多个免疫相关集合较高。癌旁仍可能含炎症或其他肝脏病变。</p></div></section><section class="section"><h2>实际检出的差异证据</h2><div class="table-scroll"><table><thead><tr><th>基因</th><th class="numeric">log2FC</th><th class="numeric">BH padj</th><th>较高组别</th></tr></thead><tbody>${wanted.map((g) => `<tr><td><button data-jump-gene="${esc(g.gene)}" class="table-row-button">${esc(g.gene)}</button></td><td class="numeric">${sign(g.lfc)}</td><td class="numeric">${sci(g.padj)}</td><td>${groupTag(g.lfc>0?'GroupX_higher':'GroupY_higher')}</td></tr>`).join('')}</tbody></table></div></section><section class="section"><h2>解读的边界</h2><p class="evidence-text">表达较高不等于存在突变。通路名称与某种癌症有关，也不能把样本诊断成该癌种。这次耐药相关 KEGG 通路没有达到校正后显著阈值；癌旁较高的药物代谢特征不能直接当作肿瘤耐药证据。</p><p class="evidence-text">免疫富集来自 bulk RNA-seq，会混合肿瘤、肝细胞、免疫与基质等讯号。仅靠本次富集不能确定免疫细胞比例、细胞活性或因果机制。</p><p class="evidence-text">完成表达与通路判定后，以 TCGA 样本类型独立核对，8 / 8 个样本与判定一致。这项核对是分析后的验证。</p></section>${renderSources()}`;
  }
  function renderSources() {
    const src = D.sources || D.references || [];
    if (!src.length) return '';
    return `<section class="section"><h2>方法与生物学数据来源</h2><ol class="source-list">${src.map((s) => typeof s === 'string' ? `<li>${esc(s)}</li>` : `<li><a href="${esc(s.url || s.link || '#')}" target="_blank" rel="noopener">${esc(s.title || s.citation || s.name || s.url)}</a>${s.detail?`　${esc(s.detail)}`:''}${s.doi?`　DOI ${esc(s.doi)}`:''}</li>`).join('')}</ol></section>`;
  }

  function render() {
    hideTooltip(); canvasObserver?.disconnect(); volcano = null; redraw = null;
    const step = D.steps.find((s) => s.id === ids[active]) || D.steps[active];
    $('step-count').textContent = `${active+1} / 9`;
    $('step-number').textContent = `STEP ${String(active+1).padStart(2,'0')}`;
    $('step-topic').textContent = topics[active];
    $('step-nav').innerHTML = ids.map((id,i) => `${i%3===0?`<div class="nav-phase">${['数据与设计','样本与基因','通路与判定'][i/3]}</div>`:''}<a href="#${id}" class="step-link ${i===active?'active':''}" ${i===active?'aria-current="step"':''}><span class="step-index">${String(i+1).padStart(2,'0')}</span><span>${navNames[i]}</span>${i<active?'<span class="step-done" aria-hidden="true">✓</span>':''}</a>`).join('');
    const widgets = [overview,quality,design,sampleStructure,differential,heatmap,() => enrichment('ora'),() => enrichment('gsea'),interpretation];
    $('step-content').innerHTML = active===0 ? widgets[0]() : `<header class="article-heading"><h1>${esc(step.title || navNames[active])}</h1><p class="intro">${esc(step.summary || '')}</p></header>${widgets[active]()}${renderFindings(step)}${step.note?`<p class="method-note">${esc(step.note)}</p>`:''}${figureButtons(step.figures)}`;
    document.body.classList.toggle('overview',active===0);
    $('journey-progress').style.width=`${(active+1)/9*100}%`;
    if(active===8 && D.interpretation?.genes?.length){
      const section=document.createElement('section');section.className='section';
      section.innerHTML=`<h2>区分癌基因作用与肿瘤标记</h2><p class="plot-hint">展开查看文献功能分类；下列机制来自外部研究，并非本次表达数据的功能验证。</p>${D.interpretation.genes.map((g)=>`<details class="gene-function"><summary><strong>${esc(g.gene)}</strong><span>${esc(g.classification)}</span></summary><p>${esc(g.explanation)}${g.source?.url?` <a href="${esc(g.source.url)}" target="_blank" rel="noopener">${esc(g.source.title)}${g.source.detail?`（${esc(g.source.detail)}）`:''}</a>`:''}</p></details>`).join('')}`;
      const source=document.querySelector('.source-list')?.closest('section');
      $('step-content').insertBefore(section,source || document.querySelector('.result-reading'));
    }
    $('previous').disabled = active===0; $('next').disabled = active===8;
    $('navigation-label').textContent = `${String(active+1).padStart(2,'0')} / 09`;
    $('next').innerHTML = active === 8 ? '已读完' : `下一步 <span aria-hidden="true">→</span>`;
    renderCode();
    if (active === 3) drawSamples();
    if (active === 4) initVolcano();
    if (active === 5) drawHeatmap();
    if (active === 6 || active === 7) {bindEnrichment();drawEnrichment();}
    document.title = `${navNames[active]} · TCGA–LIHC 配对分析`;
  }
  function setStep(index) {
    const next = Math.max(0,Math.min(8,index));
    if (next === active) return;
    location.hash = ids[next];
  }
  function readHash() {
    const [id,query] = location.hash.slice(1).split('?');
    active = Math.max(0,ids.indexOf(id));
    const gene = new URLSearchParams(query || '').get('gene');
    if (gene && genes.has(gene.toUpperCase())) state.gene = gene.toUpperCase();
    render();
    window.scrollTo({top:0,behavior:'instant'});
    const nav = document.querySelector('.step-link.active');
    if (nav && window.innerWidth <=700) $('step-nav').scrollLeft = Math.max(0,nav.offsetLeft - 25);
  }
  function highlight(code) {
    const re = /#[^\n]*|"(?:\\.|[^"\\])*"|'(?:\\.|[^'\\])*'|\b(?:TRUE|FALSE|NULL|NA|if|else|for|function|library|require|stopifnot)\b|\b\d+(?:\.\d+)?(?:e[-+]?\d+)?\b/g;
    let html = '', start=0;
    for (const m of code.matchAll(re)) {
      html += esc(code.slice(start,m.index));
      const cls = m[0][0]==='#'?'r-comment':/^["']/.test(m[0])?'r-string':/^\d/.test(m[0])?'r-number':'r-keyword';
      html += `<span class="${cls}">${esc(m[0])}</span>`; start = m.index+m[0].length;
    }
    return html+esc(code.slice(start));
  }
  function renderCode() {
    const step = D.steps.find((s) => s.id === ids[active]) || D.steps[active];
    currentCode = codeMode==='step' ? step.code || '# 请下载完整 R 分析代码。' : runCode;
    $('code-text').innerHTML = currentCode.split('\n').map((line,i)=>`<span class="code-line"><span class="code-line-number" aria-hidden="true">${i+1}</span><span class="code-line-content">${highlight(line)||' '}</span></span>`).join('');
    $('code-heading-title').textContent=`${String(active+1).padStart(2,'0')} · ${navNames[active]}`;
    $('code-source').textContent = codeMode === 'step' ? step.codeSource || 'TCGA_LIHC_analysis.R' : '在本地 R 环境执行';
    $('code-tab').setAttribute('aria-selected',String(codeMode === 'step'));
    $('run-tab').setAttribute('aria-selected',String(codeMode === 'run'));
    $('code-pre').scrollTop = 0;
  }
  function setCodeOpen(open) {
    $('code-panel').classList.toggle('mobile-open',open);
    $('code-panel').setAttribute('aria-hidden',String(!open));
    $('code-toggle').setAttribute('aria-expanded',String(open));
    $('code-backdrop').hidden=!open;
    document.body.classList.toggle('code-is-open',open);
    (open?$('code-close'):$('code-toggle')).focus();
  }
  function toast(message) { $('toast').textContent=message; $('toast').classList.add('visible');clearTimeout(toast.timer);toast.timer=setTimeout(() => $('toast').classList.remove('visible'),2200); }
  async function copyCode() {
    try {
      if(navigator.clipboard && window.isSecureContext) {
        await navigator.clipboard.writeText(currentCode);
      } else {
        const ta=document.createElement('textarea');ta.value=currentCode;ta.style.position='fixed';ta.style.opacity='0';document.body.append(ta);ta.select();
        let copied=false;try{copied=document.execCommand('copy');}finally{ta.remove();}
        if(!copied)throw Error('Clipboard unavailable');
      }
      toast(codeMode==='step'?'已复制本步代码':'已复制运行说明');
    } catch {toast('复制未完成，请选择代码后复制。');}
  }
  function tooltip(event,html) {
    const el=$('chart-tooltip');el.innerHTML=html;el.hidden=false;
    const box=el.getBoundingClientRect();let x=event.clientX+13,y=event.clientY+12;
    if(x+box.width>window.innerWidth-12)x=event.clientX-box.width-12;
    if(y+box.height>window.innerHeight-12)y=event.clientY-box.height-12;
    el.style.left=`${Math.max(8,x)}px`;el.style.top=`${Math.max(8,y)}px`;
  }
  function hideTooltip() { $('chart-tooltip').hidden=true; }
  function openFigure(png,title,pdf) {
    $('lightbox-title').textContent=title;$('lightbox-img').src=png;$('lightbox-img').alt=title;
    $('lightbox-pdf').href=pdf || png.replace(/\.png$/i,'.pdf');
    $('lightbox').showModal();hideTooltip();
  }

  function drawSamples() {
    const pca = D.pca || {points:[],variance:[]};
    if (state.sampleTab === 'correlation') {
      $('sample-plot-title').textContent='样本表达相关性';
      drawMatrix('sample-plot',D.correlation,false);
      $('sample-readout').textContent='移到格子上查看两个样本的 Pearson r。相关性以全部 VST 基因计算。';
      return;
    }
    $('sample-plot-title').textContent='PCA，变异最高的 500 个基因';
    const points=pca.points || [],w=660,h=380,p={l:62,r:30,t:28,b:55};
    if (!points.length) {$('sample-plot').innerHTML='<p class="empty-state">没有 PCA 座标数据。</p>';return;}
    const xs=points.map((v)=>Number(v.pc1)),ys=points.map((v)=>Number(v.pc2));
    const xmin=Math.min(...xs)-10,xmax=Math.max(...xs)+10,ymin=Math.min(...ys)-9,ymax=Math.max(...ys)+9;
    const xx=(v)=>p.l+(v-xmin)/(xmax-xmin)*(w-p.l-p.r),yy=(v)=>h-p.b-(v-ymin)/(ymax-ymin)*(h-p.t-p.b);
    let svg=`<svg viewBox="0 0 ${w} ${h}" role="img" aria-label="PCA 配对样本图">`;
    for(let i=0;i<5;i++){const xv=xmin+(xmax-xmin)*i/4,yv=ymin+(ymax-ymin)*i/4;svg+=`<line x1="${xx(xv)}" y1="${p.t}" x2="${xx(xv)}" y2="${h-p.b}" stroke="#e9eef2"/><line x1="${p.l}" y1="${yy(yv)}" x2="${w-p.r}" y2="${yy(yv)}" stroke="#e9eef2"/><text x="${xx(xv)}" y="${h-p.b+20}" text-anchor="middle" fill="#81909b" font-size="10">${Math.round(xv)}</text><text x="${p.l-10}" y="${yy(yv)+3}" text-anchor="end" fill="#81909b" font-size="10">${Math.round(yv)}</text>`;}
    for(const patient of ['P1','P2','P3','P4']){const pair=points.filter((v)=>v.patient===patient);if(pair.length===2)svg+=`<line x1="${xx(pair[0].pc1)}" y1="${yy(pair[0].pc2)}" x2="${xx(pair[1].pc1)}" y2="${yy(pair[1].pc2)}" stroke="${state.patient===patient?'#263f55':'#c2cdd6'}" stroke-width="${state.patient===patient?2:1}" opacity="${state.patient&&state.patient!==patient?.2:1}"/>`;}
    points.forEach((v,i)=>{const opacity=state.patient&&state.patient!==v.patient?.18:1;const px=xx(v.pc1),py=yy(v.pc2);const dx=v.group==='GroupX'?8: (v.patient==='P2'?9:-8);const dy=v.patient==='P4'?-10:v.patient==='P3'?16:-10;svg+=`<g data-pca="${i}" tabindex="0" role="button" aria-label="${esc(v.sample)}, ${esc(v.patient)}, ${esc(v.group)}" style="cursor:pointer" opacity="${opacity}"><circle cx="${px}" cy="${py}" r="${state.patient===v.patient?7:5.5}" fill="${colors[v.group]}" stroke="#fff" stroke-width="2"/><text x="${px+dx}" y="${py+dy}" text-anchor="${dx<0?'end':'start'}" fill="${colors[v.group]}" font-size="10">${esc(v.sample)}</text></g>`;});
    svg+=`<text x="${(w+p.l-p.r)/2}" y="${h-10}" text-anchor="middle" fill="#596f80" font-size="11">PC1 (${fmt(pca.variance[0],1)}%)</text><text transform="translate(16 ${(h-p.b+p.t)/2}) rotate(-90)" text-anchor="middle" fill="#596f80" font-size="11">PC2 (${fmt(pca.variance[1],1)}%)</text></svg>`;
    $('sample-plot').innerHTML=svg;
    $('sample-plot').querySelectorAll('[data-pca]').forEach((el)=>{const v=points[Number(el.dataset.pca)];const pick=()=>{state.patient=state.patient===v.patient?null:v.patient;drawSamples();updatePatientButtons();$('sample-readout').innerHTML=`<strong>${esc(v.sample)}</strong>　${esc(v.patient)} / ${esc(v.group)}　PC1 ${fmt(v.pc1)}，PC2 ${fmt(v.pc2)}`;};el.addEventListener('click',pick);el.addEventListener('keydown',(e)=>{if(e.key==='Enter'||e.key===' '){e.preventDefault();pick();}});el.addEventListener('pointermove',(e)=>tooltip(e,`<strong>${esc(v.sample)}</strong>${esc(v.patient)}　${esc(v.group)}<br>PC1 ${fmt(v.pc1)}<br>PC2 ${fmt(v.pc2)}`));el.addEventListener('pointerleave',hideTooltip);});
  }
  function updatePatientButtons() {document.querySelectorAll('[data-patient]').forEach((b)=>b.setAttribute('aria-pressed',String((b.dataset.patient||null)===state.patient)));}
  function drawMatrix(container,matrix,isHeat) {
    if (!matrix?.samples?.length) {$(container).innerHTML='<p class="empty-state">没有矩阵数据。</p>';return;}
    const cols=matrix.samples;let colIndices=cols.map((_,i)=>i);
    if(isHeat)colIndices.sort((a,b)=>state.heatOrder==='group'?groupOf(cols[a]).localeCompare(groupOf(cols[b]))||patientOf(cols[a]).localeCompare(patientOf(cols[b])):patientOf(cols[a]).localeCompare(patientOf(cols[b]))||groupOf(cols[a]).localeCompare(groupOf(cols[b])));
    const rows=isHeat?matrix.genes:matrix.samples;
    const cellW=isHeat?60:52,cellH=isHeat?18:42,left=isHeat?105:70,top=isHeat?86:72;
    const w=left+cellW*cols.length+35,h=top+cellH*rows.length+65;
    let svg=`<svg viewBox="0 0 ${w} ${h}" style="min-width:${isHeat?610:530}px" role="img" aria-label="${isHeat?'40 基因 VST Z-score 热图':'VST 样本 Pearson 相关性矩阵'}">`;
    const patientPal={P1:'#427c76',P2:'#ad7897',P3:'#ba932a',P4:'#5b9bc2'};
    colIndices.forEach((ci,j)=>{const id=cols[ci],x=left+j*cellW;svg+=`<rect x="${x}" y="${top-30}" width="${cellW-1}" height="8" fill="${patientPal[patientOf(id)]||'#8a9daa'}"/><rect x="${x}" y="${top-18}" width="${cellW-1}" height="9" fill="${colors[groupOf(id)]}"/><text x="${x+cellW/2}" y="${top-39}" text-anchor="middle" fill="#4a6478" font-size="10">${esc(id)}</text>`;});
    svg+=`<text x="${left-9}" y="${top-23}" text-anchor="end" font-size="9" fill="#6a8090">患者</text><text x="${left-9}" y="${top-10}" text-anchor="end" font-size="9" fill="#6a8090">组别</text>`;
    rows.forEach((row,i)=>{svg+=`<text x="${left-9}" y="${top+i*cellH+cellH/2+3.5}" text-anchor="end" fill="#476074" font-size="${isHeat?10:11}">${esc(row)}</text>`;colIndices.forEach((ci,j)=>{const value=Number(matrix.values[i][ci]);const fill=isHeat?heatColor(value):mix('#f6f9fb','#0072b2',Math.max(0,Math.min(1,(value-.6)/.4)));const x=left+j*cellW,y=top+i*cellH;svg+=`<g data-matrix-i="${i}" data-matrix-j="${ci}"><rect x="${x}" y="${y}" width="${cellW-1}" height="${cellH-1}" fill="${fill}"/>${!isHeat?`<text x="${x+cellW/2}" y="${y+cellH/2+3}" text-anchor="middle" font-size="9" fill="${value>.88?'#fff':'#4e6475'}" pointer-events="none">${value.toFixed(3)}</text>`:''}</g>`;});});
    svg+='</svg>';$(container).innerHTML=svg;
    $(container).querySelectorAll('[data-matrix-i]').forEach((el)=>{const i=Number(el.dataset.matrixI),j=Number(el.dataset.matrixJ),v=Number(matrix.values[i][j]);el.addEventListener('pointermove',(e)=>tooltip(e,isHeat?`<strong>${esc(rows[i])}</strong>${esc(cols[j])}　${esc(patientOf(cols[j]))} / ${groupOf(cols[j])}<br>逐基因 Z-score：${v.toFixed(3)}`:`<strong>${esc(rows[i])} × ${esc(cols[j])}</strong>Pearson r：${v.toFixed(4)}`));el.addEventListener('pointerleave',hideTooltip);});
  }
  function mix(a,b,t){const c=(s)=>[1,3,5].map((i)=>parseInt(s.slice(i,i+2),16));const aa=c(a),bb=c(b);return `rgb(${aa.map((v,i)=>Math.round(v+(bb[i]-v)*t)).join(',')})`;}
  function heatColor(v){return mix('#f8fafb',v>=0?'#d55e00':'#0072b2',Math.min(1,Math.abs(v)/1.8));}
  function drawHeatmap(){drawMatrix('heatmap-plot',D.heatmap,true);}

  function initVolcano(){
    const canvas=$('volcano');
    volcano={canvas,points:[],hover:null};
    $('fdr-slider').addEventListener('input',(e)=>{state.fdr=Number(e.target.value);drawVolcano();updateGeneInspector();});
    $('lfc-slider').addEventListener('input',(e)=>{state.lfc=Number(e.target.value);drawVolcano();updateGeneInspector();});
    $('reset-thresholds').addEventListener('click',()=>{state.fdr=.05;state.lfc=1;$('fdr-slider').value=.05;$('lfc-slider').value=1;drawVolcano();updateGeneInspector();});
    $('gene-search').addEventListener('input',geneSuggestions);
    $('gene-search').addEventListener('keydown',(e)=>{if(e.key==='Enter'){const query=e.target.value.trim().toUpperCase();const match=genes.get(query)||de.find((g)=>g.gene.toUpperCase().startsWith(query));if(match)selectGene(match.gene);else toast('受检基因中找不到这个名称。');}});
    canvas.addEventListener('pointermove',(e)=>{const rect=canvas.getBoundingClientRect();const x=e.clientX-rect.left,y=e.clientY-rect.top;let best=null,dist=70;for(const p of volcano.points){const d=(p.x-x)**2+(p.y-y)**2;if(d<dist){dist=d;best=p;}}volcano.hover=best;canvas.style.cursor=best?'pointer':'default';if(best)tooltip(e,`<strong>${esc(best.g.gene)}</strong>log2FC：${sign(best.g.lfc)}<br>BH padj：${sci(best.g.padj)}<br><span class="tt-muted">点击固定查看</span>`);else hideTooltip();});
    canvas.addEventListener('pointerleave',hideTooltip);
    canvas.addEventListener('click',()=>{if(volcano.hover)selectGene(volcano.hover.g.gene);});
    redraw=drawVolcano;drawVolcano();updateGeneInspector();
    if(window.ResizeObserver){canvasObserver=new ResizeObserver(()=>drawVolcano());canvasObserver.observe(canvas);}
  }
  function geneSuggestions(){
    const query=$('gene-search').value.trim().toUpperCase(),el=$('gene-suggestions');
    if(!query){el.hidden=true;return;}
    const hits=de.filter((g)=>g.gene.toUpperCase().includes(query)).sort((a,b)=>(a.gene.toUpperCase()===query?-1:b.gene.toUpperCase()===query?1:0)||a.gene.localeCompare(b.gene)).slice(0,7);
    el.innerHTML=hits.length?hits.map((g)=>`<button data-gene="${esc(g.gene)}">${esc(g.gene)} <span class="${g.lfc>0?'x-text':'y-text'}">${sign(g.lfc)}</span></button>`).join(''):'<div style="padding:9px 12px;font-size:11px">受检基因中没有相符名称</div>';el.hidden=false;
  }
  function selectGene(gene){
    state.gene=gene.toUpperCase();$('gene-search').value=gene;$('gene-suggestions').hidden=true;
    history.replaceState(null,'',`#differential?gene=${encodeURIComponent(gene)}`);drawVolcano();updateGeneInspector();hideTooltip();
  }
  function updateGeneInspector(){
    const g=genes.get(state.gene);if(!g){$('gene-inspector').innerHTML='<p>搜索或点击基因，固定查看数值。</p>';return;}
    const hit=finite(g.padj)&&g.padj<state.fdr&&Math.abs(g.lfc)>state.lfc;
    $('gene-inspector').innerHTML=`<div><strong class="gene-name">${esc(g.gene)}</strong><div class="gene-subtitle">${hit?'符合当前探索阈值':'未达当前探索阈值'}</div></div><div><small>GroupX / GroupY 的 log2FC</small><span class="number ${g.lfc>0?'x-text':'y-text'}">${sign(g.lfc)}</span><small>${g.lfc>0?'GroupX 较高':'GroupY 较高'}</small></div><div><small>BH 校正后 P 值</small><span class="number">${sci(g.padj)}</span><small>baseMean ${fmt(g.baseMean,1)}</small></div>`;
  }
  function drawVolcano(){
    if(!volcano||!$('volcano'))return;
    const {canvas}=volcano;const width=canvas.clientWidth,height=canvas.clientHeight,dpr=window.devicePixelRatio||1;
    canvas.width=Math.round(width*dpr);canvas.height=Math.round(height*dpr);const c=canvas.getContext('2d');c.scale(dpr,dpr);
    const pad={l:52,r:19,t:22,b:47},plotW=width-pad.l-pad.r,plotH=height-pad.t-pad.b;
    const valid=de.filter((g)=>finite(g.padj)&&finite(g.lfc));
    const xmin=Math.floor(Math.min(-5,...valid.map((g)=>g.lfc)))-.6,xmax=Math.ceil(Math.max(5,...valid.map((g)=>g.lfc)))+.6;
    const maxY=Math.ceil(Math.max(10,...valid.map((g)=>-Math.log10(Math.max(Number(g.padj),1e-300))))/10)*10;
    const xx=(x)=>pad.l+(x-xmin)/(xmax-xmin)*plotW,yy=(y)=>height-pad.b-y/maxY*plotH;
    c.fillStyle='#fff';c.fillRect(0,0,width,height);c.font='10px Arial';c.textAlign='right';c.lineWidth=1;
    for(let i=0;i<=4;i++){const y=maxY*i/4;c.strokeStyle='#eaf0f4';c.beginPath();c.moveTo(pad.l,yy(y));c.lineTo(width-pad.r,yy(y));c.stroke();c.fillStyle='#83929e';c.fillText(String(Math.round(y)),pad.l-8,yy(y)+3);}
    c.textAlign='center';for(let x=Math.ceil(xmin/5)*5;x<xmax;x+=5){c.strokeStyle='#eef2f5';c.beginPath();c.moveTo(xx(x),pad.t);c.lineTo(xx(x),height-pad.b);c.stroke();c.fillStyle='#83929e';c.fillText(String(x),xx(x),height-pad.b+18);}
    c.setLineDash([4,4]);c.strokeStyle='#8c9ca8';for(const x of [-state.lfc,state.lfc]){c.beginPath();c.moveTo(xx(x),pad.t);c.lineTo(xx(x),height-pad.b);c.stroke();}c.beginPath();c.moveTo(pad.l,yy(-Math.log10(state.fdr)));c.lineTo(width-pad.r,yy(-Math.log10(state.fdr)));c.stroke();c.setLineDash([]);
    let nx=0,ny=0;volcano.points=[];const selected=genes.get(state.gene);
    for(const g of valid){const significant=g.padj<state.fdr&&Math.abs(g.lfc)>state.lfc;let color='#b9c5ce';if(significant){color=g.lfc>0?colors.GroupX:colors.GroupY;g.lfc>0?nx++:ny++;}const x=xx(g.lfc),y=yy(-Math.log10(Math.max(Number(g.padj),1e-300)));c.fillStyle=color;c.globalAlpha=significant?.72:.40;c.beginPath();c.arc(x,y,significant?1.8:1.4,0,Math.PI*2);c.fill();volcano.points.push({x,y,g});}
    c.globalAlpha=1;
    if(selected&&finite(selected.padj)){const x=xx(selected.lfc),y=yy(-Math.log10(Math.max(selected.padj,1e-300)));c.strokeStyle=colors[selected.lfc>0?'GroupX':'GroupY'];c.lineWidth=1.6;c.beginPath();c.arc(x,y,6,0,Math.PI*2);c.stroke();c.font='bold 11px Arial';c.textAlign=x>width*.75?'right':'left';const tx=x+(x>width*.75?-10:10),ty=Math.max(pad.t+12,y-8);c.fillStyle='#fff';const tw=c.measureText(selected.gene).width;c.fillRect(tx+(x>width*.75?-tw:0)-3,ty-11,tw+6,15);c.fillStyle='#1e3d56';c.fillText(selected.gene,tx,ty);}
    c.fillStyle='#5f7485';c.font='10px Arial';c.textAlign='center';c.fillText('log2 fold change (positive = GroupX higher)',pad.l+plotW/2,height-10);c.save();c.translate(13,pad.t+plotH/2);c.rotate(-Math.PI/2);c.fillText('-log10(BH-adjusted P)',0,0);c.restore();
    $('count-x').textContent=fmt(nx,0);$('count-y').textContent=fmt(ny,0);$('count-other').textContent=fmt(de.length-nx-ny,0);$('fdr-value').textContent=state.fdr.toFixed(3);$('lfc-value').textContent=state.lfc.toFixed(1);$('threshold-state').textContent=Math.abs(state.fdr-.05)<1e-8&&Math.abs(state.lfc-1)<1e-8?'正式阈值':'探索阈值';
  }

  function bindEnrichment(){
    $('term-search').addEventListener('input',(e)=>{enrichmentState.query=e.target.value;enrichmentState.page=0;drawEnrichment();});
    $('sig-only').addEventListener('change',(e)=>{enrichmentState.sig=e.target.checked;enrichmentState.page=0;drawEnrichment();});
    $('table-prev').addEventListener('click',()=>{enrichmentState.page--;drawEnrichment();});
    $('table-next').addEventListener('click',()=>{enrichmentState.page++;drawEnrichment();});
  }
  function drawEnrichment(){
    const st=enrichmentState,isGsea=st.kind==='gsea',all=(isGsea?D.gsea:D.ora)||[];
    const directionRows=all.filter((r)=>r.database===st.db&&(st.direction==='all'||r.direction===st.direction));
    let rows=directionRows.filter((r)=>(!st.sig||(finite(r.fdr)&&r.fdr<.05))&&(!st.query||`${r.id} ${r.description}`.toLowerCase().includes(st.query.toLowerCase()))).sort((a,b)=>(finite(a.fdr)?a.fdr:Infinity)-(finite(b.fdr)?b.fdr:Infinity));
    const nPages=Math.max(1,Math.ceil(rows.length/30));st.page=Math.max(0,Math.min(st.page,nPages-1));
    const shown=rows.slice(st.page*30,st.page*30+30);
    $('enrichment-status').textContent=`${st.db==='GO_BP'?'GO BP':'KEGG'}　${fmt(rows.length,0)} / ${fmt(directionRows.length,0)} 个集合`;
    $('table-page').textContent=`${st.page+1} / ${nPages}`;$('table-prev').disabled=st.page===0;$('table-next').disabled=st.page===nPages-1;
    $('enrichment-table').innerHTML=`<thead><tr><th>通路／生物过程</th><th>方向</th>${isGsea?'<th class="numeric">NES</th>':'<th class="numeric">命中基因</th>'}<th class="numeric">BH ${isGsea?'FDR':'padj'}</th></tr></thead><tbody>${shown.length?shown.map((r)=>`<tr><td class="term-name">${esc(r.description)}<small>${esc(r.id)}</small></td><td data-label="方向">${groupTag(r.direction)}</td><td class="numeric" data-label="${isGsea?'NES':'命中基因'}">${isGsea?sign(r.nes):fmt(r.count,0)}</td><td class="numeric" data-label="BH ${isGsea?'FDR':'padj'}">${fdrTag(r.fdr)}</td></tr>`).join(''):'<tr><td colspan="4"><div class="empty-state">没有符合筛选条件的集合。可修改搜索词或取消「只看 BH &lt; 0.05」。</div></td></tr>'}</tbody>`;
    const top=directionRows.filter((r)=>finite(r.fdr)&&r.fdr<.05).sort((a,b)=>a.fdr-b.fdr).slice(0,8);
    if(!top.length)$('enrichment-chart').innerHTML='<div class="empty-state">这个数据库与方向没有通过 BH &lt; 0.05 的富集结果。</div>';
    else{
      const max=Math.max(...top.map((r)=>isGsea?Math.abs(r.nes):-Math.log10(Math.max(r.fdr,1e-300))));
      $('enrichment-chart').innerHTML=`<div class="section-head"><h2>${isGsea?'排序富集':'显著富集'}的前 8 个集合</h2><small>固定 BH &lt; 0.05</small></div><div class="term-bars">${top.map((r)=>{const value=isGsea?r.nes:-Math.log10(Math.max(r.fdr,1e-300));return `<div class="term-bar" title="${esc(r.description)}&#10;${esc(r.id)}&#10;BH: ${sci(r.fdr)}${isGsea?`&#10;NES: ${sign(r.nes)}`:''}"><div class="term-bar-label">${esc(r.description)}<small>${esc(r.id)}</small></div><div class="term-bar-track"><i style="width:${Math.max(1,Math.abs(value)/max*100)}%;background:${colors[r.direction]}"></i></div><span class="bar-value">${isGsea?sign(value):value.toFixed(1)}</span></div>`;}).join('')}</div><div class="bar-axis-note">${isGsea?'条长表示 |NES|，正负号与色彩表示方向':'条长表示 −log10(BH 校正 P 值)'}</div>`;
    }
    const curves=$('gsea-curves');if(curves){
      const fs=figures.filter((f)=>f.png?.includes(`${st.db}_GSEA_curve_`)&&(st.direction==='all'||f.png.includes(st.direction)));
      curves.innerHTML=isGsea&&fs.length?`<section class="section"><h2>代表性富集曲线</h2><div class="curve-grid">${fs.map((f)=>`<div><button class="image-preview" data-figure="${esc(f.png)}" data-title="${esc(f.title)}" data-pdf="${esc(f.pdf||'')}"><img src="${esc(f.png)}" alt="${esc(f.title)}" loading="lazy"></button><div class="curve-caption">${esc(f.title)}，点图放大并查看正式输出。</div></div>`).join('')}</div></section>`:'';
    }
  }

  document.addEventListener('click',(e)=>{
    if(e.target.closest('[data-open-code]')){setCodeOpen(true);return;}
    const figure=e.target.closest('[data-figure]');if(figure){openFigure(figure.dataset.figure,figure.dataset.title,figure.dataset.pdf);return;}
    const gene=e.target.closest('[data-gene]');if(gene){selectGene(gene.dataset.gene);return;}
    const jump=e.target.closest('[data-jump-gene]');if(jump){location.hash=`differential?gene=${encodeURIComponent(jump.dataset.jumpGene)}`;return;}
    const patient=e.target.closest('[data-patient]');if(patient){state.patient=patient.dataset.patient||null;drawSamples();updatePatientButtons();return;}
    const sampleTab=e.target.closest('[data-sample-tab]');if(sampleTab){state.sampleTab=sampleTab.dataset.sampleTab;document.querySelectorAll('[data-sample-tab]').forEach((b)=>b.setAttribute('aria-pressed',String(b===sampleTab)));drawSamples();return;}
    const heatOrder=e.target.closest('[data-heat-order]');if(heatOrder){state.heatOrder=heatOrder.dataset.heatOrder;document.querySelectorAll('[data-heat-order]').forEach((b)=>b.setAttribute('aria-pressed',String(b===heatOrder)));drawHeatmap();return;}
    const db=e.target.closest('[data-db]');if(db){enrichmentState.db=db.dataset.db;enrichmentState.page=0;document.querySelectorAll('[data-db]').forEach((b)=>b.setAttribute('aria-pressed',String(b===db)));drawEnrichment();return;}
    const direction=e.target.closest('[data-direction]');if(direction){enrichmentState.direction=direction.dataset.direction;enrichmentState.page=0;document.querySelectorAll('[data-direction]').forEach((b)=>b.setAttribute('aria-pressed',String(b===direction)));drawEnrichment();return;}
    if(!$('gene-search')?.contains(e.target)&&!$('gene-suggestions')?.contains(e.target)&&$('gene-suggestions'))$('gene-suggestions').hidden=true;
  });
  $('previous').addEventListener('click',()=>setStep(active-1));$('next').addEventListener('click',()=>setStep(active+1));
  $('code-tab').addEventListener('click',()=>{codeMode='step';renderCode();});$('run-tab').addEventListener('click',()=>{codeMode='run';renderCode();});$('copy-code').addEventListener('click',copyCode);
  $('code-toggle').addEventListener('click',()=>setCodeOpen($('code-panel').getAttribute('aria-hidden')==='true'));
  $('code-close').addEventListener('click',()=>setCodeOpen(false));
  $('code-backdrop').addEventListener('click',()=>setCodeOpen(false));
  $('close-lightbox').addEventListener('click',()=>$('lightbox').close());
  $('lightbox').addEventListener('click',(e)=>{if(e.target===$('lightbox'))$('lightbox').close();});
  document.addEventListener('keydown',(e)=>{if(e.key==='Tab'&&$('code-panel').classList.contains('mobile-open')){const targets=[...$('code-panel').querySelectorAll('button,a[href],pre[tabindex]')];const first=targets[0],last=targets[targets.length-1];if(e.shiftKey&&document.activeElement===first){e.preventDefault();last.focus();}else if(!e.shiftKey&&document.activeElement===last){e.preventDefault();first.focus();}return;}if(e.key==='Escape'){if($('code-panel').classList.contains('mobile-open'))setCodeOpen(false);hideTooltip();return;}if($('lightbox').open||$('code-panel').classList.contains('mobile-open')||/INPUT|TEXTAREA|SELECT|BUTTON/.test(e.target.tagName)||e.ctrlKey||e.metaKey||e.altKey)return;if(e.key==='ArrowRight'){e.preventDefault();setStep(active+1);}if(e.key==='ArrowLeft'){e.preventDefault();setStep(active-1);}});
  window.addEventListener('hashchange',readHash);window.addEventListener('resize',()=>{clearTimeout(window.tcgaResizeTimer);window.tcgaResizeTimer=setTimeout(()=>redraw?.(),100);});
  window.addEventListener('scroll',hideTooltip,{passive:true});
  if(window.innerWidth<=1080)$('code-toggle').setAttribute('aria-expanded','false');
  readHash();
})();
