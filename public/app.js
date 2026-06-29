const form = document.querySelector("#analysis-form");
const fileInput = document.querySelector("#vcf-file");
const dropZone = document.querySelector("#drop-zone");
const fileLabel = document.querySelector("#file-label");
const sampleSelect = document.querySelector("#sample-name");
const errorBox = document.querySelector("#form-error");
const analyzeButton = document.querySelector("#analyze-button");
const resultsSection = document.querySelector("#results");
const summaryGrid = document.querySelector("#summary-grid");
const findingsList = document.querySelector("#findings-list");
const scoresList = document.querySelector("#scores-list");
const warningsList = document.querySelector("#warnings-list");
const findingCount = document.querySelector("#finding-count");
const exportButton = document.querySelector("#export-button");
let selectedFile = null;
let latestReport = null;

const element = (tag, className, text) => {
  const node = document.createElement(tag);
  if (className) node.className = className;
  if (text !== undefined) node.textContent = text;
  return node;
};

const formatBytes = (bytes) => {
  if (bytes < 1024) return `${bytes} B`;
  if (bytes < 1024 ** 2) return `${(bytes / 1024).toFixed(1)} KB`;
  return `${(bytes / 1024 ** 2).toFixed(1)} MB`;
};

const showError = (message) => {
  errorBox.textContent = message;
  errorBox.hidden = false;
};

const clearError = () => {
  errorBox.hidden = true;
  errorBox.textContent = "";
};

async function readSampleNames(file) {
  const preview = await file.slice(0, 1024 * 1024).text();
  const header = preview.split(/\r?\n/).find((line) => line.startsWith("#CHROM"));
  return header ? header.split("\t").slice(9).filter(Boolean) : [];
}

async function selectFile(file) {
  clearError();
  if (!file) return;
  if (!/\.vcf(\.gz)?$/i.test(file.name)) {
    showError("Choose a file ending in .vcf or .vcf.gz.");
    return;
  }
  if (file.size > 25 * 1024 * 1024) {
    showError("The selected file exceeds the 25 MB upload limit.");
    return;
  }
  selectedFile = file;
  fileLabel.textContent = file.name;
  dropZone.querySelector("#file-help").textContent = `${formatBytes(file.size)} · click to replace`;
  dropZone.classList.add("has-file");
  sampleSelect.replaceChildren(new Option("First available sample", ""));
  sampleSelect.disabled = true;

  if (!file.name.toLowerCase().endsWith(".gz")) {
    const names = await readSampleNames(file);
    for (const name of names) sampleSelect.add(new Option(name, name));
    sampleSelect.disabled = names.length === 0;
  }
}

dropZone.addEventListener("click", () => fileInput.click());
dropZone.addEventListener("keydown", (event) => {
  if (event.key === "Enter" || event.key === " ") {
    event.preventDefault();
    fileInput.click();
  }
});
fileInput.addEventListener("change", () => selectFile(fileInput.files[0]));
for (const eventName of ["dragenter", "dragover"]) {
  dropZone.addEventListener(eventName, (event) => {
    event.preventDefault();
    dropZone.classList.add("dragging");
  });
}
for (const eventName of ["dragleave", "drop"]) {
  dropZone.addEventListener(eventName, (event) => {
    event.preventDefault();
    dropZone.classList.remove("dragging");
  });
}
dropZone.addEventListener("drop", (event) => selectFile(event.dataTransfer.files[0]));

function renderSummary(report) {
  const values = [
    [report.summary.recordsProcessed, "Records processed"],
    [report.summary.exactAlternateMatches, "Exact matches"],
    [report.summary.highPriority, "Strong evidence"],
    [report.summary.reviewRequired, "Calls to review"]
  ];
  summaryGrid.replaceChildren(...values.map(([value, label]) => {
    const item = element("div", "summary-item");
    item.append(element("strong", "", String(value)), element("span", "", label));
    return item;
  }));
}

function tag(text, variant = "") {
  return element("span", `tag ${variant}`.trim(), text);
}

function detail(label, value) {
  const wrapper = element("div");
  wrapper.append(element("span", "", label), element("strong", "", value ?? "Not available"));
  return wrapper;
}

function renderFindings(report) {
  findingCount.textContent = `${report.findings.length} result${report.findings.length === 1 ? "" : "s"}`;
  if (!report.findings.length) {
    findingsList.replaceChildren(element(
      "div",
      "empty-state",
      "No exact alternate-allele matches were found in the educational database. This does not mean no clinically relevant variants are present."
    ));
    return;
  }

  findingsList.replaceChildren(...report.findings.map((finding) => {
    const card = element("article", "finding");
    const top = element("div", "finding-top");
    const heading = element("div");
    heading.append(element("h4", "", finding.disease), element("p", "finding-description", finding.description));
    const tags = element("div", "tags");
    tags.append(
      tag(finding.classification, finding.priority),
      tag(finding.callQuality.status === "pass" ? "Call passed checks" : "Review call", finding.callQuality.status)
    );
    top.append(heading, tags);
    const metadata = element("div", "finding-details");
    metadata.append(
      detail("Locus", `${finding.locus} ${finding.ref}>${finding.alt}`),
      detail("Genotype", `${finding.genotype} · ${finding.zygosity}`),
      detail("Gene", finding.gene),
      detail("Molecular consequence", finding.geneImpact),
      detail(
        "Evidence",
        finding.evidence.status !== "unverified"
          ? `${finding.evidence.sourceName ?? "Source matched"} · ${finding.evidence.reviewStars}/4 review stars · ${finding.evidence.accession}`
          : "Unverified project annotation"
      ),
      detail("Site / genotype quality", `${finding.callQuality.issues.length ? finding.callQuality.issues.join(", ") : "No configured flags"}`),
      detail("Database odds ratio", String(finding.oddsRatio)),
      detail("Source", `VCF line ${finding.sourceLine}`)
    );
    if (finding.activeModifiers.length) {
      const modifierText = finding.activeModifiers
        .map((item) => `${item.factor} (database multiplier ${item.databaseMultiplier})`).join(", ");
      metadata.append(detail("Active annotations", modifierText));
    }
    if (finding.ageAnnotation) {
      const ageText = finding.ageAnnotation.status === "reference-only"
        ? `${Math.round(finding.ageAnnotation.databasePenetrance * 100)}% at age ${finding.ageAnnotation.age} — reference annotation only`
        : `No database point before age ${finding.ageAnnotation.firstReferenceAge}`;
      metadata.append(detail("Age annotation", ageText));
    }
    card.append(top, metadata);
    return card;
  }));
}

function renderScores(report) {
  if (!report.relativeScores.length) {
    scoresList.replaceChildren(element("div", "warning", "Risk scoring is disabled; this report shows classifications and call quality only."));
    return;
  }
  const maximum = Math.max(...report.relativeScores.map((item) => Math.abs(item.logOddsScore)), 1);
  scoresList.replaceChildren(...report.relativeScores.map((score) => {
    const item = element("div", "score");
    item.append(
      element("strong", "", score.disease),
      element("span", "", `Relative OR multiplier ${score.relativeOddsMultiplier} · ${score.contributors} contributor(s)`)
    );
    const bar = element("div", "score-bar");
    const fill = element("i");
    fill.style.width = `${Math.max(4, Math.abs(score.logOddsScore) / maximum * 100)}%`;
    bar.append(fill);
    item.append(bar);
    return item;
  }));
}

function renderWarnings(report) {
  if (!report.warnings.length) {
    warningsList.replaceChildren(element("div", "warning", "No parser warnings."));
    return;
  }
  warningsList.replaceChildren(...report.warnings.slice(0, 20).map((warning) => {
    const item = element("div", "warning");
    item.append(
      element("strong", "", warning.code.replaceAll("_", " ")),
      element("span", "", `${warning.line ? `Line ${warning.line}: ` : ""}${warning.message}`)
    );
    return item;
  }));
}

function renderReport(report) {
  latestReport = report;
  renderSummary(report);
  document.querySelector("#report-notice").textContent = report.disclaimer;
  renderFindings(report);
  renderScores(report);
  renderWarnings(report);
  resultsSection.hidden = false;
  resultsSection.scrollIntoView({ behavior: "smooth", block: "start" });
}

form.addEventListener("submit", async (event) => {
  event.preventDefault();
  clearError();
  if (!selectedFile) {
    showError("Choose a VCF file before starting the analysis.");
    dropZone.focus();
    return;
  }

  const age = document.querySelector("#age").value;
  if (age && (!Number.isInteger(Number(age)) || Number(age) < 1 || Number(age) > 120)) {
    showError("Age must be a whole number between 1 and 120.");
    return;
  }

  const payload = new FormData();
  payload.append("vcf", selectedFile);
  if (age) payload.append("age", age);
  if (sampleSelect.value) payload.append("sampleName", sampleSelect.value);
  const genomeBuild = document.querySelector("#genome-build").value;
  if (genomeBuild) payload.append("genomeBuild", genomeBuild);
  payload.append("factors", "[]");

  analyzeButton.disabled = true;
  analyzeButton.querySelector("span").textContent = "Analyzing…";
  try {
    const response = await fetch("/api/v1/analyze", { method: "POST", body: payload });
    const data = await response.json().catch(() => null);
    if (!response.ok) throw new Error(data?.error?.message ?? "The analysis request failed.");
    renderReport(data);
  } catch (error) {
    showError(error.message || "The analysis could not be completed.");
  } finally {
    analyzeButton.disabled = false;
    analyzeButton.querySelector("span").textContent = "Review my file";
  }
});

exportButton.addEventListener("click", () => {
  if (!latestReport) return;
  const blob = new Blob([JSON.stringify(latestReport, null, 2)], { type: "application/json" });
  const link = document.createElement("a");
  link.href = URL.createObjectURL(blob);
  link.download = `helix-report-${new Date().toISOString().slice(0, 10)}.json`;
  link.click();
  URL.revokeObjectURL(link.href);
});
