package uk.ac.ebi.spot.goci.ui.controller;

import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.http.HttpEntity;
import org.springframework.http.MediaType;
import org.springframework.stereotype.Controller;
import org.springframework.ui.Model;
import org.springframework.web.bind.annotation.PathVariable;
import org.springframework.web.bind.annotation.RequestMapping;
import org.springframework.web.bind.annotation.RequestMethod;
import org.springframework.web.bind.annotation.RequestParam;
import org.springframework.web.bind.annotation.ResponseBody;
import org.springframework.web.client.RestTemplate;
import org.springframework.web.util.UriComponentsBuilder;
import uk.ac.ebi.spot.goci.model.EfoTraitResult;
import uk.ac.ebi.spot.goci.ui.SearchConfiguration;

/**
 * Created by xinhe on 19/04/2017.
 */
@Controller
public class EfoTraitController {


    private SearchConfiguration searchConfiguration;

    @Autowired
    public EfoTraitController(SearchConfiguration searchConfiguration) {
        this.searchConfiguration = searchConfiguration;
    }


//    @RequestMapping(value = "/efotrait2/{efoId}", produces = MediaType.TEXT_HTML_VALUE)
//    public String search2(Model model,
//                         @PathVariable(required = false) String efoId,
//                         @RequestParam(required = false) String filter) {
//        EfoTraitResult efoTraitResult = new EfoTraitResult();
//        efoTraitResult.setQuery(efoId);
//        efoTraitResult.setFilter(filter);
//        efoTraitResult.setEfoId(efoId);
//        model.addAttribute("result", efoTraitResult);
//        return "efotrait-page-rest";
//    }

    @RequestMapping(value = "efotraits", produces = MediaType.TEXT_HTML_VALUE)
    public String search() {
        return "efotraits";
    }

    @RequestMapping(value = "efotraits/{efoId}", produces = MediaType.TEXT_HTML_VALUE)
    public String search(Model model,
                         @PathVariable(required = false) String efoId,
                         @RequestParam(required = false) String filter,
                         @RequestParam(required = false, defaultValue = "") String included,
                         @RequestParam(required = false, defaultValue = "") String checked) {
        EfoTraitResult efoTraitResult = new EfoTraitResult();
        efoTraitResult.setQuery(efoId);
        efoTraitResult.setFilter(filter);
        efoTraitResult.setEfoId(efoId);
        efoTraitResult.setIncluded(included);
        efoTraitResult.setChecked(checked);
        model.addAttribute("result", efoTraitResult);
        return "efotrait-page";
    }

//    @RequestMapping(value = "/trait/test", method = RequestMethod.GET)
//    public @ResponseBody String searchSolr(@RequestParam(required = false) String filter) {
//        EfoTraitResult efoTraitResult = new EfoTraitResult();
//        efoTraitResult.setQuery("EFO_0000400");
//        efoTraitResult.setFilter(filter);
//        efoTraitResult.setEfoId("EFO_0000400");
//
//        RestTemplate restTemplate = new RestTemplate();
//        UriComponentsBuilder builder = UriComponentsBuilder.fromHttpUrl("http://localhost:8280/gwas/api/search/efotrait")
//                .queryParam("q","shortForm:EFO_0001645")
//                .queryParam("fq","shortForm:EFO_0001645")
//                .queryParam("facet","on")
//                .queryParam("facrt.field","resourcename")
//                .queryParam("fl","shortForm");
//        HttpEntity<String> response= restTemplate.getForEntity(builder.build().encode().toUri(), String.class);
//        return response.getBody().toString();
//    }


}
