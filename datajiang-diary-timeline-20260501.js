const START_DATE = new Date(2025, 0, 1);
const END_DATE = new Date(2026, 4, 1);
const DAY_MS = 24 * 60 * 60 * 1000;

// 构建后的页面会从 window.DIARY_ENTRIES 读取真实内容；没有注入数据时使用这里的示例。
const demoEntries = [
  { date: "2025-01-01", type: "day", text: "新的一年开始，先把心放稳，把日子写清楚。" },
  { date: "2025-02-01", type: "month", text: "二月像一次短暂停靠，整理、休息，也继续观察自己。" },
  { date: "2025-05-01", type: "day", text: "今天适合慢下来，看见自己已经走了不少路。" },
  { date: "2025-10-01", type: "month", text: "十月把一些旧问题翻出来，也把答案往前推了一点。" },
  { date: "2026-01-01", type: "day", text: "重新开头，不必喊得很响，先认真过好今天。" },
  { date: "2026-05-01", type: "day", text: "这一段长时间的记录先收束成一条能回看的线。" },
];
const demoNodes = demoEntries.map((entry) => ({ date: entry.date }));
const customEntries = Array.isArray(window.DIARY_ENTRIES) ? window.DIARY_ENTRIES : demoEntries;
const customNodes = Array.isArray(window.DIARY_NODES) ? window.DIARY_NODES : demoNodes;

const typeText = {
  day: "日记",
  month: "月记",
};

const state = {
  mode: "all",
  query: "",
  focusEntryId: null,
  focusStack: [],
  panorama: false,
  panShiftX: null,
  panAnimation: null,
  layout: null,
};

const elements = {
  track: document.querySelector("#timelineTrack"),
  pointLayer: document.querySelector("#pointLayer"),
  cardLayer: document.querySelector("#cardLayer"),
  monthBands: document.querySelector("#monthBands"),
  statusLine: document.querySelector("#statusLine"),
  backButton: document.querySelector("#backButton"),
  panLeftButton: document.querySelector("#panLeftButton"),
  panRightButton: document.querySelector("#panRightButton"),
  panoramaButton: document.querySelector("#panoramaButton"),
  searchInput: document.querySelector("#searchInput"),
  buttons: [...document.querySelectorAll(".mode-button")],
};

const allEntries = buildEntries();
const allNodes = buildNodes(allEntries);
render();

elements.buttons.forEach((button) => {
  button.addEventListener("click", () => {
    cancelPanAnimation();
    state.mode = button.dataset.mode;
    clearFocus();
    state.panorama = false;
    state.panShiftX = null;
    elements.buttons.forEach((item) => item.classList.toggle("active", item === button));
    render();
  });
});

elements.searchInput.addEventListener("input", (event) => {
  cancelPanAnimation();
  state.query = event.target.value.trim().toLowerCase();
  clearFocus();
  state.panorama = false;
  state.panShiftX = null;
  render();
});

elements.pointLayer.addEventListener("click", (event) => {
  const point = event.target.closest(".point");

  if (!point) {
    return;
  }

  pushFocus(point.dataset.entryId, event.clientX);
  render();
});

elements.cardLayer.addEventListener("click", (event) => {
  const card = event.target.closest(".entry-card");

  if (!card) {
    return;
  }

  pushFocus(card.dataset.nodeId, event.clientX);
  render();
});

elements.backButton.addEventListener("click", () => {
  cancelPanAnimation();
  state.focusStack.pop();
  state.focusEntryId = state.focusStack.at(-1)?.id ?? null;
  state.panShiftX = null;
  render();
});

elements.panLeftButton.addEventListener("click", () => panTimeline(1));
elements.panRightButton.addEventListener("click", () => panTimeline(-1));
elements.panoramaButton.addEventListener("click", () => {
  cancelPanAnimation();
  clearFocus();
  state.panorama = true;
  state.panShiftX = null;
  render();
});

window.addEventListener("resize", () => {
  window.requestAnimationFrame(render);
});

function buildEntries() {
  return customEntries
    .filter((entry) => entry && entry.date && ["day", "month"].includes(entry.type))
    .map((entry) => ({
      id: `${entry.date}:${entry.type}`,
      date: entry.date,
      type: entry.type,
      text: String(entry.text ?? ""),
    }))
    .sort((left, right) => parseISO(left.date) - parseISO(right.date) || typeRank(left.type) - typeRank(right.type));
}

function buildNodes(entries) {
  const nodeMap = new Map();

  buildDateRange(getRangeStart(entries), getRangeEnd(entries)).forEach((date) => {
    nodeMap.set(date, { id: `${date}:node`, date });
  });

  customNodes.forEach((node) => {
    if (node?.date) {
      nodeMap.set(node.date, { id: `${node.date}:node`, date: node.date });
    }
  });

  entries.forEach((entry) => {
    nodeMap.set(entry.date, { id: `${entry.date}:node`, date: entry.date });
  });

  return [...nodeMap.values()].sort((left, right) => parseISO(left.date) - parseISO(right.date));
}

function buildDateRange(startDate, endDate) {
  const dates = [];
  const current = parseISO(startDate);
  const end = parseISO(endDate);

  while (current <= end) {
    dates.push(toISO(current));
    current.setDate(current.getDate() + 1);
  }

  return dates;
}

function getRangeStart(entries) {
  const dates = [
    toISO(START_DATE),
    ...customNodes.map((node) => node?.date).filter(Boolean),
    ...entries.map((entry) => entry.date),
  ];

  return dates.sort((left, right) => left.localeCompare(right))[0];
}

function getRangeEnd(entries) {
  const dates = [
    toISO(END_DATE),
    ...customNodes.map((node) => node?.date).filter(Boolean),
    ...entries.map((entry) => entry.date),
  ];

  return dates.sort((left, right) => right.localeCompare(left))[0];
}

function render() {
  const filtered = getFilteredEntries();
  const nodes = getFilteredNodes(filtered);
  state.layout = getLayout(nodes);
  const visibleNodes = getVisibleNodes(nodes, state.layout);
  const visibleEntries = getVisibleEntries(filtered, state.layout);
  drawTrack(state.layout);
  drawMonthLabels(visibleEntries, state.layout);
  drawPoints(visibleNodes, state.layout);
  drawCards(visibleEntries, state.layout);
  updateBackButton();
  updatePanButtons();
  updateStatus(filtered, nodes);
}

function getLayout(nodes) {
  const width = window.innerWidth;
  const height = window.innerHeight;
  const mobile = width < 760;
  const left = mobile ? 26 : Math.max(42, width * 0.035);
  const right = width - (mobile ? 26 : Math.max(42, width * 0.035));
  const top = mobile ? Math.max(188, height * 0.24) : Math.max(112, height * 0.16);
  const bottom = height - (mobile ? 98 : Math.max(88, height * 0.105));
  const depth = (state.panorama ? 0 : 1) + state.focusStack.length;
  const expansion = 3 ** depth;
  const baseColumnGap = mobile ? 64 : 92;
  const baseRowGap = mobile ? 34 : 30;
  const columnGap = baseColumnGap * expansion;
  const rowGap = baseRowGap * expansion;
  const rowsPerColumn = Math.max(1, Math.floor((bottom - top) / rowGap) + 1);
  const contentHeight = (rowsPerColumn - 1) * rowGap;
  const verticalSlack = Math.max(0, bottom - top - contentHeight);
  const timelineNodes = orderNodesForTimeline(nodes);
  const indexById = new Map(timelineNodes.map((node, index) => [node.id, index]));
  const contentColumns = Math.max(1, Math.ceil(timelineNodes.length / rowsPerColumn));
  const virtualColumns = 2;
  const totalColumns = contentColumns + virtualColumns;
  const contentWidth = (totalColumns - 1) * columnGap;
  const minShiftX = Math.min(0, right - (left + contentWidth));
  const maxShiftX = 0;
  const shiftY = verticalSlack / 2;
  const anchor = state.focusStack.at(-1);
  let shiftX = minShiftX;

  if (state.panShiftX !== null) {
    shiftX = clamp(state.panShiftX, minShiftX, maxShiftX);
  } else if (anchor && indexById.has(anchor.id)) {
    const anchorIndex = indexById.get(anchor.id);
    const anchorColumn = Math.floor(anchorIndex / rowsPerColumn) + 1;
    const anchorX = left + anchorColumn * columnGap;
    shiftX = clamp(anchor.x - anchorX, minShiftX, maxShiftX);
  }

  return {
    width,
    height,
    mobile,
    left,
    right,
    top,
    bottom,
    depth,
    expansion,
    rowsPerColumn,
    columnGap,
    rowGap,
    contentHeight,
    timelineNodes,
    indexById,
    contentColumns,
    virtualColumns,
    totalColumns,
    minShiftX,
    maxShiftX,
    shiftX,
    shiftY,
  };
}

function drawTrack(layout) {
  const firstColumn = Math.max(1, Math.floor((0 - layout.left - layout.shiftX) / layout.columnGap) - 1);
  const lastColumn = Math.min(
    layout.contentColumns,
    Math.ceil((layout.width - layout.left - layout.shiftX) / layout.columnGap) + 1,
  );
  const columns = [];

  for (let index = firstColumn; index <= lastColumn; index += 1) {
    const x = layout.left + layout.columnGap * index + layout.shiftX;
    const yTop = layout.top + layout.shiftY;
    const yBottom = yTop + layout.contentHeight;
    const visibleTop = clamp(yTop, layout.top, layout.bottom);
    const visibleBottom = clamp(yBottom, layout.top, layout.bottom);
    columns.push(`
      <span class="track-column" style="left:${x}px;top:${visibleTop}px;height:${Math.max(0, visibleBottom - visibleTop)}px"></span>
      <span class="track-cap start" style="left:${x}px;top:${visibleTop}px"></span>
      <span class="track-cap end" style="left:${x}px;top:${visibleBottom}px"></span>
    `);
  }

  elements.track.innerHTML = columns.join("");
}

function drawMonthLabels(entries, layout) {
  const months = entries.filter((entry) => entry.type === "month");
  elements.monthBands.innerHTML = months
    .map((entry, index) => {
      const point = positionForEntry(entry, layout);
      const side = point.column % 2 === 0 ? 1 : -1;
      const x = clamp(point.x + side * 30, 36, layout.width - 36);
      const y = clamp(point.y - 18 + (index % 2) * 12, layout.top - 28, layout.bottom + 24);
      return `<span class="month-label" style="left:${x}px;top:${y}px">${formatMonth(entry.date)}</span>`;
    })
    .join("");
}

function drawPoints(entries, layout) {
  elements.pointLayer.innerHTML = entries
    .map((node) => {
      const point = positionForNode(node, layout);
      const entriesAtNode = allEntries.filter((entry) => entry.date === node.date);
      const hasMonth = entriesAtNode.some((entry) => entry.type === "month");
      const pointType = hasMonth ? "month" : entriesAtNode.length > 0 ? "day" : "empty";
      const active = state.focusEntryId === node.id ? " active" : "";
      return `<button class="point ${pointType}${active}" data-entry-id="${node.id}" style="left:${point.x}px;top:${point.y}px" type="button" aria-label="${formatDate(node.date)} 节点"></button>`;
    })
    .join("");
}

function drawCards(entries, layout) {
  const orderedEntries = orderCardsForLayer(entries);
  elements.cardLayer.innerHTML = orderedEntries
    .map((entry, index) => {
      const point = cardPositionFor(entry, index, layout);
      const active = state.focusEntryId === entry.id || state.focusEntryId === `${entry.date}:node` ? " active" : "";
      return `
        <article class="entry-card ${entry.type}${active}" data-entry-id="${entry.id}" data-node-id="${entry.date}:node" style="left:${point.x}px;top:${point.y}px">
          <div class="entry-meta">
            <time datetime="${entry.date}">${formatDate(entry.date)}</time>
            <span class="entry-type">${typeText[entry.type]}</span>
          </div>
          <p>${escapeHTML(entry.text)}</p>
        </article>
      `;
    })
    .join("");
}

function getFilteredEntries() {
  return allEntries.filter((entry) => {
    const modeMatch = state.mode === "all" || entry.type === state.mode;
    const queryMatch =
      !state.query ||
      entry.date.includes(state.query) ||
      entry.text.toLowerCase().includes(state.query) ||
      typeText[entry.type].includes(state.query);
    return modeMatch && queryMatch;
  });
}

function getFilteredNodes(entries) {
  const entriesByDate = new Map();

  entries.forEach((entry) => {
    if (!entriesByDate.has(entry.date)) {
      entriesByDate.set(entry.date, []);
    }

    entriesByDate.get(entry.date).push(entry);
  });

  return allNodes.filter((node) => {
    const entriesAtNode = entriesByDate.get(node.date) ?? [];
    const modeMatch = state.mode === "all" || entriesAtNode.length > 0;
    const queryMatch =
      !state.query ||
      node.date.includes(state.query) ||
      entriesAtNode.some((entry) => entry.text.toLowerCase().includes(state.query) || typeText[entry.type].includes(state.query));
    return modeMatch && queryMatch;
  });
}

function getVisibleNodes(nodes, layout) {
  const marginX = layout.mobile ? 190 : 250;
  const marginY = layout.mobile ? 90 : 110;

  return nodes.filter((node) => {
    const point = positionForNode(node, layout);
    return (
      point.x >= -marginX &&
      point.x <= layout.width + marginX &&
      point.y >= layout.top - marginY &&
      point.y <= layout.bottom + marginY
    );
  });
}

function getVisibleEntries(entries, layout) {
  const marginX = layout.mobile ? 190 : 250;
  const marginY = layout.mobile ? 90 : 110;

  return entries.filter((entry) => {
    const point = positionForEntry(entry, layout);
    return (
      point.x >= -marginX &&
      point.x <= layout.width + marginX &&
      point.y >= layout.top - marginY &&
      point.y <= layout.bottom + marginY
    );
  });
}

function orderCardsForLayer(entries) {
  const lowerCards = entries.filter((entry) => entry.type !== "month");
  const monthCards = entries.filter((entry) => entry.type === "month");
  return [...lowerCards, ...monthCards];
}

function orderEntriesForTimeline(entries) {
  return [...entries].sort((left, right) => parseISO(left.date) - parseISO(right.date) || typeRank(left.type) - typeRank(right.type));
}

function orderNodesForTimeline(nodes) {
  return [...nodes].sort((left, right) => parseISO(left.date) - parseISO(right.date));
}

function positionForNode(node, layout) {
  const index = layout.indexById.get(node.id) ?? 0;
  const column = Math.floor(index / layout.rowsPerColumn) + 1;
  const row = index % layout.rowsPerColumn;
  const x = layout.left + column * layout.columnGap + layout.shiftX;
  const y = layout.top + row * layout.rowGap + layout.shiftY;

  return { x, y, column, row, index };
}

function positionForEntry(entry, layout) {
  return positionForNode({ id: `${entry.date}:node`, date: entry.date }, layout);
}

function cardPositionFor(entry, index, layout) {
  const base = positionForEntry(entry, layout);
  const typeNudge = entry.type === "month" ? -10 : 10;
  let x = base.x;
  const y = base.y + typeNudge;

  if (layout.depth === 0) {
    x += ((index % 3) - 1) * 10;
  }

  return { x, y };
}

function panTimeline(direction) {
  const filtered = getFilteredEntries();
  const layout = state.layout ?? getLayout(getFilteredNodes(filtered));
  const step = Math.max(layout.columnGap, (layout.right - layout.left) * 0.72);
  const targetShiftX = clamp(layout.shiftX + direction * step, layout.minShiftX, layout.maxShiftX);
  animatePanTo(targetShiftX);
}

function animatePanTo(targetShiftX) {
  const filtered = getFilteredEntries();
  const layout = state.layout ?? getLayout(getFilteredNodes(filtered));
  const startShiftX = layout.shiftX;
  const distance = Math.abs(targetShiftX - startShiftX);

  if (distance < 1) {
    return;
  }

  if (state.panAnimation) {
    window.cancelAnimationFrame(state.panAnimation.frame);
  }

  const duration = clamp(260 + distance * 0.18, 320, 620);
  const startedAt = performance.now();

  const tick = (now) => {
    const progress = Math.min(1, (now - startedAt) / duration);
    const eased = easeInOutCubic(progress);
    state.panShiftX = startShiftX + (targetShiftX - startShiftX) * eased;
    render();

    if (progress < 1) {
      state.panAnimation = { frame: window.requestAnimationFrame(tick) };
    } else {
      state.panShiftX = targetShiftX;
      state.panAnimation = null;
      render();
    }
  };

  state.panAnimation = { frame: window.requestAnimationFrame(tick) };
}

function easeInOutCubic(value) {
  return value < 0.5 ? 4 * value ** 3 : 1 - (-2 * value + 2) ** 3 / 2;
}

function updateBackButton() {
  const hasHistory = state.focusStack.length > 0;
  elements.backButton.classList.toggle("visible", hasHistory);
  elements.backButton.disabled = !hasHistory;
}

function updatePanButtons() {
  const layout = state.layout;
  const canPanLeft = layout.shiftX < layout.maxShiftX - 1;
  const canPanRight = layout.shiftX > layout.minShiftX + 1;
  elements.panLeftButton.disabled = !canPanLeft;
  elements.panRightButton.disabled = !canPanRight;
}

function updateStatus(entries, nodes) {
  const modeName = state.mode === "all" ? "全部" : typeText[state.mode];
  const query = state.query ? `，搜索：${state.query}` : "";
  const focus = allEntries.find((entry) => entry.id === state.focusEntryId);
  const depth = state.layout?.depth ?? 0;
  const focusText = focus
    ? `，第 ${depth} 层展开 ${formatDate(focus.date)}`
    : state.panorama
      ? "，第 0 层全景"
      : "，默认第 1 层显示最近记录";
  elements.statusLine.value = `${modeName}${query}，当前显示 ${nodes.length} 个节点，${entries.length} 条记录${focusText}。`;
}

function pushFocus(entryId, x = window.innerWidth / 2) {
  if (!entryId) {
    return;
  }

  cancelPanAnimation();
  state.panShiftX = null;
  state.focusStack.push({ id: entryId, x });
  state.focusEntryId = entryId;
}

function clearFocus() {
  state.focusEntryId = null;
  state.focusStack = [];
}

function cancelPanAnimation() {
  if (!state.panAnimation) {
    return;
  }

  window.cancelAnimationFrame(state.panAnimation.frame);
  state.panAnimation = null;
}

function toISO(date) {
  const year = date.getFullYear();
  const month = String(date.getMonth() + 1).padStart(2, "0");
  const day = String(date.getDate()).padStart(2, "0");
  return `${year}-${month}-${day}`;
}

function parseISO(value) {
  const [year, month, day] = value.split("-").map(Number);
  return new Date(year, month - 1, day);
}

function formatDate(value) {
  const date = parseISO(value);
  return `${date.getFullYear()}.${date.getMonth() + 1}.${date.getDate()}`;
}

function formatMonth(value) {
  const date = parseISO(value);
  return `${date.getFullYear()}.${String(date.getMonth() + 1).padStart(2, "0")}`;
}

function typeRank(type) {
  return type === "month" ? 0 : 1;
}

function clamp(value, min, max) {
  return Math.min(max, Math.max(min, value));
}

function escapeHTML(value) {
  return value
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;")
    .replaceAll('"', "&quot;")
    .replaceAll("'", "&#039;");
}
