(() => {
  const select=document.getElementById('report-page');
  const previous=document.getElementById('report-previous'),next=document.getElementById('report-next');
  let page=1;
  select.innerHTML=Array.from({length:17},(_,i)=>`<option value="${i+1}">${i+1}</option>`).join('');
  function show(n){
    page=n;select.value=String(n);previous.disabled=n===1;next.disabled=n===17;
    document.getElementById('report-img').src=`assets/report/page-${n}.webp`;
    document.getElementById('report-img').alt=`扩展报告第 ${n} 页`;
    document.getElementById('report-image-link').href=`assets/report/page-${n}.webp`;
    document.getElementById('report-status').textContent=`第 ${n} 页 / 17 · 点击页面可放大阅读`;
  }
  select.onchange=()=>show(Number(select.value));previous.onclick=()=>show(page-1);next.onclick=()=>show(page+1);
})();
